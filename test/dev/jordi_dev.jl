
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

# This respects the input triangulation, without adding back all the ghosts
function MyFESpace(
  trian::GridapDistributed.DistributedTriangulation,reffe;split_own_and_ghost=false,constraint=nothing,kwargs...
)
  spaces = map(local_views(trian)) do t
    FESpace(t,reffe;kwargs...)
  end
  gids = GridapDistributed.generate_gids(trian,spaces)
  vector_type = GridapDistributed._find_vector_type(spaces,gids;split_own_and_ghost=split_own_and_ghost)
  space = GridapDistributed.DistributedSingleFieldFESpace(spaces,gids,trian,vector_type)
  return GridapDistributed._add_distributed_constraint(space,reffe,constraint)
end

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
trian = Triangulation(single_ghost, agg_cutgeo, ACTIVE)
writevtk(agg_model,"data/agg_model");
writevtk(Triangulation(no_ghost, agg_cutgeo, ACTIVE),"data/agg_trian");
map(ranks,local_views(trian)) do r, trian
  writevtk(trian,"data/agg_trian_$(r)");
end

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V = MyFESpace(trian,reffe)
V = MyFESpace(trian,reffe;dirichlet_tags=["boundary"])

# (TMP) This gives wrong output when root is not owned, but does not affect below.
lcell_to_lroot_bis = _local_aggregates(lcell_to_root,PRange(cell_indices))
@assert lcell_to_lroot_bis == lcell_to_lroot

bgcell_to_bgroot = lcell_to_lroot
Vagg = AgFEMSpace(V,bgcell_to_bgroot)
