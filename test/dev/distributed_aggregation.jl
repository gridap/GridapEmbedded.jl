module DistributedAggregationTest

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays

using GridapEmbedded: aggregate
using GridapEmbedded.Distributed: _local_aggregates

using Test
using BenchmarkTools

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

function _global_aggregates(cell_to_lcellin,gids::PRange)
  map(_global_aggregates,cell_to_lcellin,local_to_global(gids))
end

function _global_aggregates(cell_to_lcellin,lcell_to_gcell)
  T = eltype(cell_to_lcellin)
  map(cell_to_lcellin) do lcin
    iszero(lcin) ? lcin : T(lcell_to_gcell[ lcin ])
  end
end

function run_old_distributed_aggregation(ranks,
                                         parts,
                                         ncells_x_dir,
                                         problem)

  bgmodel, geo = problem(ranks, parts, ncells_x_dir, 1)
  cutgeo = cut(bgmodel, geo)
  strategy = AggregateCutCellsByThreshold(1.0)
  bgmodel,_,lcell_to_root = aggregate(strategy,cutgeo)
  bgmodel,lcell_to_root
end

function run_new_distributed_aggregation(ranks,
                                         parts,
                                         ncells_x_dir,
                                         nghost_layers,
                                         problem)

  bgmodel, geo = problem(ranks, parts, ncells_x_dir, nghost_layers)
  cutgeo = cut(bgmodel, geo)

  gids = get_cell_gids(bgmodel)
  cell_indices = partition(gids)

  strategy = AggregateCutCellsByThreshold(1.0)
  lcell_to_lroot, lcell_to_root, lcell_to_value =
    map(local_views(cutgeo),cell_indices) do cutgeo,cell_indices
      lid_to_gid = local_to_global(cell_indices)
      aggregate(strategy,cutgeo,geo,lid_to_gid,IN)
    end |> tuple_of_arrays

  lcell_to_owner = map(copy∘local_to_owner,cell_indices)
  lcell_to_owner = map(lcell_to_owner,lcell_to_lroot) do lcell_to_owner,lcell_to_lroot
    for i in eachindex(lcell_to_owner)
      if !iszero(lcell_to_lroot[i])
        lcell_to_owner[i] = lcell_to_owner[lcell_to_lroot[i]]
      end
    end
    lcell_to_owner
  end

  lcell_to_root,_ =
    find_optimal_roots!(lcell_to_root,lcell_to_value,lcell_to_owner,cell_indices);

  bgmodel,lcell_to_root
end

function run_benchmark_test(distribute,
                            parts,
                            ncells_x_dir,
                            nghost_layers,
                            problem)

  ranks = distribute(LinearIndices((prod(parts),)))
  @btime obgmodel,olcell_to_root = run_old_distributed_aggregation(
                                    $ranks,parts,ncells_x_dir,problem)
  @btime nbgmodel,nlcell_to_root = run_new_distributed_aggregation(
                                    $ranks,parts,ncells_x_dir,nghost_layers,problem)

  # ogids = get_cell_gids(obgmodel)
  # olcell_to_root = _global_aggregates(olcell_to_root,ogids)
  # oocell_to_root = map(olcell_to_root,own_to_local(ogids)) do agg,o_to_l
  #   map(Reindex(agg),o_to_l)
  # end

  # ngids = get_cell_gids(nbgmodel)
  # nocell_to_root = map(nlcell_to_root,own_to_local(ngids)) do agg,o_to_l
  #   map(Reindex(agg),o_to_l)
  # end

  # writevtk(
  #   Triangulation(obgmodel), "data/kettlebell_aggregates_old", 
  #   celldata = ["aggregate" => oocell_to_root],
  # );

  # writevtk(
  #   Triangulation(nbgmodel), "data/kettlebell_aggregates_new", 
  #   celldata = ["aggregate" => nocell_to_root],
  # );

  # # Note: The criterion to choose the root is slightly 
  # # different between the old and new aggregation. For
  # # this reason, the following test will generally fail.
  # # In contrast to the old aggregation, the new aggregation
  # # chooses the root cell with lowest GID among the ones
  # # that are reached at the same num. of iterations and
  # # have the same centroid distance.
  # map(oocell_to_root,nocell_to_root) do oocell_to_root,nocell_to_root
  #   @test all(oocell_to_root .== nocell_to_root)
  # end

end

function run_new_distributed_agfem(distribute,
                                   parts,
                                   ncells_x_dir,
                                   nghost_layers,
                                   problem)

  ranks = distribute(LinearIndices((prod(parts),)))
  bgmodel, geo = problem(ranks, parts, ncells_x_dir, nghost_layers)
  cutgeo = cut(bgmodel, geo)

  gids = get_cell_gids(bgmodel)
  cell_indices = partition(gids)

  strategy = AggregateCutCellsByThreshold(1.0)
  lcell_to_lroot, lcell_to_root, lcell_to_value =
    map(local_views(cutgeo),cell_indices) do cutgeo,cell_indices
      lid_to_gid = local_to_global(cell_indices)
      aggregate(strategy,cutgeo,geo,lid_to_gid,IN)
    end |> tuple_of_arrays

  lcell_to_owner = map(copy∘local_to_owner,cell_indices)
  lcell_to_owner = map(lcell_to_owner,lcell_to_lroot) do lcell_to_owner,lcell_to_lroot
    for i in eachindex(lcell_to_owner)
      if !iszero(lcell_to_lroot[i])
        lcell_to_owner[i] = lcell_to_owner[lcell_to_lroot[i]]
      end
    end
    lcell_to_owner
  end

  lcell_to_inconsistent_root = map(copy,lcell_to_root)

  lcell_to_root, lcell_to_owner, lcell_to_value, agg_cell_indices =
    find_optimal_roots!(lcell_to_root,lcell_to_value,lcell_to_owner,cell_indices);

  # Output for verification of lcell_to_root map

  ocell_to_root = map(lcell_to_root,own_to_local(gids)) do agg,o_to_l
    map(Reindex(agg),o_to_l)
  end

  writevtk(EmbeddedBoundary(cutgeo),"data/bnd");
  writevtk(
    Triangulation(bgmodel), "data/kettlebell_aggregates", 
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
      Triangulation(bgmodel), "data/kettlebell_aggregates_$(r)", 
      celldata = [ "roots"              => lcell_to_root,
                   "inconsistent roots" => lcell_to_iroot,
                   "owners"             => lcell_to_owner ],
    );
  end

  # "Creating" aggregate-conforming cell partition

  # RMK: Current distributed AgFEM is constructed on 
  # the local portion and assumes the root cells of  
  # all ghost cells are in the local portion.

  # For now, extract only owned cells (no_ghost).
  function active_aggregate_conforming_trian(model,cutgeo,agg_cell_indices)
    trians = map( local_views(cutgeo),
                  local_views(agg_cell_indices) ) do cutgeo, agg_cell_indices
      Triangulation(no_ghost,agg_cell_indices,cutgeo,ACTIVE)
    end
    GridapDistributed.DistributedTriangulation(trians,model)
  end

  aggtrian = active_aggregate_conforming_trian(bgmodel,cutgeo,agg_cell_indices)

  writevtk(aggtrian,"data/kettlebell_aggtrian");

  order = 1
  reffe = ReferenceFE(lagrangian,Float64,order)

  # (TMP) This gives wrong output when root is not owned, but does not affect below.
  lcell_to_lroot = _local_aggregates(lcell_to_root,gids)

  spaces = map( local_views(aggtrian),
                local_views(lcell_to_lroot),
                local_views(gids) ) do aggtrian, lcell_to_lroot, gids
    space = TestFESpace(aggtrian,reffe)
    AgFEMSpace(space,lcell_to_lroot,space,local_to_global(gids))
  end

  map(ranks,local_views(spaces)) do r,space
    aggtrian = get_triangulation(space)
    writevtk(
      aggtrian,
      "data/kettlebell_aggtrian_$(r)",
      celldata=[
        "part" => fill(r,num_cells(aggtrian)),
        # "dof_ids" => string.(get_cell_dof_ids(space))
      ],
    );
  end;

end

distribute = PartitionedArrays.DebugArray
parts = (3,1)
ncells_x_dir = 9
nghost_layers = 2
problem = symmetric_kettlebell

run_benchmark_test(distribute,
                   parts,
                   ncells_x_dir,
                   nghost_layers,
                   problem)

end