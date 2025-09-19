
using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using MPI

using Gridap.Arrays, Gridap.FESpaces
using Gridap.Geometry, Gridap.CellData

using GridapEmbedded.AgFEM, GridapEmbedded.Distributed
using GridapEmbedded: aggregate
using GridapEmbedded.Distributed: _local_aggregates, DistributedEmbeddedDiscretization

using GridapDistributed: DistributedDiscreteModel

function asymmetric_kettlebell(ranks, parts, nc, ng)
  L = 1
  p0 = Point(0.0,0.0)
  pmin = Point(-L/4,-L/8)
  pmax = Point(L/4,3*L/8)

  R = 0.2
  x = VectorValue(1.0,0.0)
  h = L/(2*nc)

  geo1 = disk(R,x0=p0)
  geo2 = quadrilateral(;x0=Point(-L/4+h/4,5*h/4),
                        d1=VectorValue(6*h+h/2,0.0),
                        d2=VectorValue(0.0,4*h+h/2))
  geo3 = quadrilateral(;x0=Point(-L/4+3*h/4,7*h/4),
                        d1=VectorValue(6*h-h/2,0.0),
                        d2=VectorValue(0.0,4*h-h/2))
  geo = union(geo1,intersect(geo2,!geo3))
  bgmodel = CartesianDiscreteModel(ranks,parts,pmin,pmax,(nc,nc);ghost=(ng,ng))
  return bgmodel, geo
end

function symmetric_kettlebell(ranks, parts, nc, ng)
  L = 1
  p0 = Point(0.0,0.0)
  pmin = Point(-L/4,-L/8)
  pmax = Point(L/4,3*L/8)

  R = 0.2
  h = L/(2*nc)

  geo1 = disk(R,x0=p0)
  geo2 = quadrilateral(;x0=Point(-L/4+h/4,5*h/4),
                        d1=VectorValue(8*h+h/2,0.0),
                        d2=VectorValue(0.0,5*h+h/2))
  geo3 = quadrilateral(;x0=Point(-L/4+3*h/4,7*h/4),
                        d1=VectorValue(8*h-h/2,0.0),
                        d2=VectorValue(0.0,5*h-h/2))
  geo = union(geo1,intersect(geo2,!geo3))
  bgmodel = CartesianDiscreteModel(ranks,parts,pmin,pmax,(nc,nc);ghost=(ng,ng))
  return bgmodel, geo
end

function exchange_impl!(vector_partition,cache)
  buffer_snd = map(vector_partition,cache) do values, cache
    local_indices_snd = cache.local_indices_snd
    for (p,lid) in enumerate(local_indices_snd.data)
      cache.buffer_snd.data[p] = values[lid]
    end
    cache.buffer_snd
  end
  neighbors_snd, neighbors_rcv, buffer_rcv = map(cache) do cache
    cache.neighbors_snd, cache.neighbors_rcv, cache.buffer_rcv
  end |> tuple_of_arrays
  graph = ExchangeGraph(neighbors_snd,neighbors_rcv)
  t = exchange!(buffer_rcv,buffer_snd,graph)
  return t
end

function find_optimal_roots!(lcell_to_root,lcell_to_value,lcell_to_owner,cell_indices)
  # Bring all root candidates to the owner of the cut cell
  roots_cache = PartitionedArrays.p_vector_cache(lcell_to_root,cell_indices)
  t1 = exchange_impl!(lcell_to_root,roots_cache)
  values_cache = PartitionedArrays.p_vector_cache(lcell_to_value,cell_indices)
  t2 = exchange_impl!(lcell_to_value,values_cache)
  owners_cache = PartitionedArrays.p_vector_cache(lcell_to_owner,cell_indices)
  t3 = exchange_impl!(lcell_to_owner,owners_cache)
  wait(t1)
  wait(t2)
  wait(t3)

  # Select the optimal root for each local cut cell
  map(
    lcell_to_root,lcell_to_value,lcell_to_owner,roots_cache,values_cache,owners_cache
  ) do lcell_to_root, lcell_to_value, lcell_to_owner, roots_cache, values_cache, owners_cache
    lids_rcv = roots_cache.local_indices_rcv
    values_rcv = values_cache.buffer_rcv
    roots_rcv = roots_cache.buffer_rcv
    owners_rcv = owners_cache.buffer_rcv
    for k in eachindex(roots_cache.neighbors_rcv)
      for (lcell,root,value,owner) in zip(lids_rcv[k], roots_rcv[k], values_rcv[k], owners_rcv[k])
        if lcell_to_value[lcell] > value # Take the minimum
          lcell_to_root[lcell] = root # Global ID
          lcell_to_value[lcell] = value
          lcell_to_owner[lcell] = owner
        end
      end
    end
  end

  # Scatter the optimal roots and values
  # Technically, we do not need to scatter the values, it can be removed after 
  # we are done debugging
  t1 = consistent!(PVector(lcell_to_owner,cell_indices,owners_cache))
  t2 = consistent!(PVector(lcell_to_root,cell_indices,roots_cache))
  t3 = consistent!(PVector(lcell_to_value,cell_indices,values_cache))
  
  wait(t1)

  agg_cell_indices = map(cell_indices,lcell_to_owner) do cell_indices, lcell_to_owner
    LocalIndices(global_length(cell_indices),part_id(cell_indices),local_to_global(cell_indices),lcell_to_owner)
  end

  wait(t2)
  wait(t3)

  return lcell_to_root, lcell_to_owner, lcell_to_value, agg_cell_indices
end

np = (2,1)
ranks = collect(1:prod(np))
bgmodel, geo = asymmetric_kettlebell(ranks, np, 8, 4)
cutgeo = cut(bgmodel, geo)

cell_indices = partition(get_cell_gids(bgmodel))

strategy = AggregateCutCellsByThreshold(1.0)
lcell_to_lroot, lcell_to_root, lcell_to_value = map(local_views(cutgeo),cell_indices) do cutgeo, cell_indices
  lid_to_gid = local_to_global(cell_indices)
  aggregate(strategy,cutgeo,geo,lid_to_gid,IN)
end |> tuple_of_arrays

lcell_to_owner = map(cell_indices,lcell_to_lroot) do cell_indices,lcell_to_lroot
  lcell_to_owner = copy(local_to_owner(cell_indices))
  for (lcell,lroot) in enumerate(lcell_to_lroot)
    if !iszero(lroot)
      lcell_to_owner[lcell] = lcell_to_owner[lroot]
    end
  end
  lcell_to_owner
end

lcell_to_inconsistent_root = map(copy,lcell_to_root)

lcell_to_root, lcell_to_owner, lcell_to_value, agg_cell_indices =
  find_optimal_roots!(lcell_to_root,lcell_to_value,lcell_to_owner,cell_indices);

# Output for verification of lcell_to_root map

ocell_to_root = map(getindex,lcell_to_root,map(own_to_local,cell_indices))
writevtk(EmbeddedBoundary(cutgeo),"data/boundary");
writevtk(
  Triangulation(bgmodel), "data/aggregates", 
  celldata = ["aggregate" => ocell_to_root],
);
map(ranks,
    local_views(bgmodel),
    lcell_to_root,
    lcell_to_inconsistent_root,
    lcell_to_owner) do r,bgmodel,
                      lcell_to_root,
                      lcell_to_iroot,
                      lcell_to_owner
  writevtk(
    Triangulation(bgmodel), "data/aggregates_$(r)", 
    celldata = [ "roots"              => lcell_to_root,
                  "inconsistent roots" => lcell_to_iroot,
                  "owners"             => lcell_to_owner ],
  );
end

agg_model = DistributedDiscreteModel(local_views(bgmodel), PRange(agg_cell_indices))
agg_cutgeo = DistributedEmbeddedDiscretization(local_views(cutgeo), agg_model)

# TODO: Remove extra layers of ghosts when building the triangulation
trian = Triangulation(agg_cutgeo, ACTIVE)
writevtk(agg_model,"data/agg_model");
writevtk(Triangulation(no_ghost, agg_cutgeo, ACTIVE),"data/agg_trian");

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V = TestFESpace(trian,reffe)

# (TMP) This gives wrong output when root is not owned, but does not affect below.
lcell_to_lroot_bis = _local_aggregates(lcell_to_root,PRange(cell_indices))
@assert lcell_to_lroot_bis == lcell_to_lroot

bgcell_to_bgroot = lcell_to_lroot
aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs, aggdof_gids, expanded_dof_gids = AgFEMSpace(V,bgcell_to_bgroot)


############################################################################################
############################################################################################
############################################################################################

f = V
bgcell_to_bgroot = lcell_to_lroot
shfns_g = get_fe_basis(f)
dofs_g = get_fe_dof_basis(f)
spaces = local_views(f)
cell_gids = GridapDistributed.generate_cell_gids(trian)

trian = get_triangulation(f)
bgcell_to_gcell = map(local_views(trian)) do trian
  glue = get_glue(trian,Val(num_cell_dims(trian)))
  glue.mface_to_tface
end

# The following is the same as in serial, it's just some reindexing onto the triangulation
nlfdofs, cell_to_root, cell_to_fdofs, cell_to_coeffs, cell_to_proj, cell_to_gcell = map(
  spaces, bgcell_to_bgroot, local_views(shfns_g), local_views(dofs_g), bgcell_to_gcell
) do space, bgcell_to_bgroot, shfns_g, dofs_g, bgcell_to_gcell
  trian = get_triangulation(space)

  glue = get_glue(trian,Val(num_cell_dims(trian)))
  cell_to_bgcell = glue.tface_to_mface
  bgcell_to_cell = glue.mface_to_tface
  cell_to_bgroot = view(bgcell_to_bgroot,cell_to_bgcell)
  cell_to_root = collect(Int32,lazy_map(Reindex(bgcell_to_cell),cell_to_bgroot))
  cell_to_gcell = collect(Int32,lazy_map(Reindex(bgcell_to_gcell),cell_to_bgcell))

  cell_phys_shfns_g = get_array(change_domain(shfns_g,PhysicalDomain()))
  cell_phys_root_shfns_g = lazy_map(Reindex(cell_phys_shfns_g),cell_to_root)
  root_shfns_g = GenericCellField(cell_phys_root_shfns_g,trian,PhysicalDomain())

  # Compute data needed to compute the constraints
  dofs_f = get_fe_dof_basis(space)
  shfns_f = get_fe_basis(space)
  cell_to_coeffs = dofs_f(root_shfns_g)
  cell_to_proj = dofs_g(shfns_f)
  cell_to_fdofs = get_cell_dof_ids(space)
  nlfdofs = num_free_dofs(space)

  return nlfdofs, cell_to_root, cell_to_fdofs, cell_to_coeffs, cell_to_proj, cell_to_gcell
end |> tuple_of_arrays

fdof_is_agg, fdof_to_cell, fdof_to_ldof = map(
  nlfdofs, cell_to_root, cell_to_fdofs, cell_to_gcell
) do nlfdofs, cell_to_root, cell_to_fdofs, cell_to_gcell
  fdof_is_agg, fdof_to_cell, fdof_to_ldof = AgFEM._allocate_fdof_to_data(nlfdofs)
  AgFEM._fill_fdof_to_data!(fdof_is_agg,fdof_to_cell,fdof_to_ldof,cell_to_root,cell_to_fdofs,cell_to_gcell)
  return fdof_is_agg, fdof_to_cell, fdof_to_ldof
end |> tuple_of_arrays

mdof_gids, aggdof_gids, mdof_to_fdof, aggdof_to_fdof, fdof_to_posneg = GridapEmbedded.Distributed.generate_aggregated_gids(
  cell_gids, cell_to_fdofs, fdof_is_agg
)

# Create aggdof to nldofs mapping
ptrs, aggdof_to_nldofs = map(
  aggdof_to_fdof,fdof_to_cell,cell_to_root,cell_to_fdofs,partition(aggdof_gids)
) do aggdof_to_fdof,fdof_to_cell,cell_to_root,cell_to_fdofs,aggdof_ids
  ptrs = GridapEmbedded.Distributed.get_aggdof_ptrs(
    aggdof_to_fdof,fdof_to_cell,cell_to_root,cell_to_fdofs,own_to_local(aggdof_ids)
  )
  aggdof_to_nldofs = view(ptrs,2:length(ptrs))
  return ptrs, aggdof_to_nldofs
end |> tuple_of_arrays
consistent!(PVector(aggdof_to_nldofs, partition(aggdof_gids))) |> wait

aggdof_to_dofs, aggdof_to_coeffs = map(
  ptrs,aggdof_to_fdof,fdof_to_cell,fdof_to_ldof,cell_to_root,cell_to_fdofs,cell_to_coeffs,cell_to_proj,partition(aggdof_gids)
) do ptrs,aggdof_to_fdof,fdof_to_cell,fdof_to_ldof,cell_to_root,cell_to_fdofs,cell_to_coeffs,cell_to_proj,aggdof_ids
  aggdof_to_dofs_data, aggdof_to_coeffs_data = AgFEM._allocate_aggdof_to_data(ptrs,cell_to_coeffs)
  GridapEmbedded.Distributed.aggdof_to_dofs!(
    aggdof_to_dofs_data,ptrs,aggdof_to_fdof,fdof_to_cell,cell_to_root,cell_to_fdofs,own_to_local(aggdof_ids)
  )
  GridapEmbedded.Distributed.aggdof_to_coeffs!(
    aggdof_to_coeffs_data,ptrs,aggdof_to_fdof,fdof_to_cell,fdof_to_ldof,cell_to_coeffs,cell_to_proj,own_to_local(aggdof_ids)
  )
  return JaggedArray(aggdof_to_dofs_data,ptrs), JaggedArray(aggdof_to_coeffs_data,ptrs)
end |> tuple_of_arrays

map(GridapEmbedded.Distributed.to_global_dofs!, aggdof_to_dofs, partition(mdof_gids), fdof_to_posneg)
t1 = consistent!(PVector(aggdof_to_dofs, partition(aggdof_gids)))
t2 = consistent!(PVector(aggdof_to_coeffs, partition(aggdof_gids)))
wait(t1)
expanded_dof_gids = map(GridapEmbedded.Distributed.to_local_dofs!,aggdof_to_dofs, partition(aggdof_gids), partition(mdof_gids), mdof_to_fdof, nlfdofs)
wait(t2)

aggdof_to_dofs, aggdof_to_coeffs = map(aggdof_to_dofs, aggdof_to_coeffs) do dofs, coeffs
  Table(dofs.data,dofs.ptrs), Table(coeffs.data,coeffs.ptrs)
end |> tuple_of_arrays
