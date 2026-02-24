module DistributedAggregation

using Gridap
using GridapEmbedded
using GridapDistributed

using PartitionedArrays
const PArrays = PartitionedArrays

using GridapEmbedded: aggregate
using GridapEmbedded.Distributed: _local_aggregates

import GridapEmbedded.LevelSetCutters: disk, sphere, popcorn, Leaf

using MPI

using Test
using BenchmarkTools

function disk(ranks, parts, nc, ng)
  geo = disk(1.0)
  pmin = Point(-1.1,-1.1)
  pmax = Point(1.1,1.1)
  bgmodel = CartesianDiscreteModel(ranks,parts,pmin,pmax,(nc,nc);ghost=(ng,ng))
  return bgmodel, geo
end

function sphere(ranks, parts, nc, ng)
  geo = sphere(1.0)
  pmin = Point(-1.1,-1.1,-1.1)
  pmax = Point(1.1,1.1,1.1)
  bgmodel = CartesianDiscreteModel(ranks,parts,pmin,pmax,(nc,nc,nc);ghost=(ng,ng,ng))
  return bgmodel, geo
end

function flower(ranks, parts, nc, ng;
    x₀=Point(0.0,0.0), R₀=0.6, m=0.6, ω=5.0)
  name="flower"
  function flowerfun(x)
    _flower(x,x₀,R₀,m,ω)
  end
  tree = Leaf((flowerfun,name,nothing))
  geo = AnalyticalGeometry(tree)
  pmin = Point(-1.1,-1.1)
  pmax = Point(1.1,1.1)
  bgmodel = CartesianDiscreteModel(ranks,parts,pmin,pmax,(nc,nc);ghost=(ng,ng))
  return bgmodel, geo
end

@inline function _flower(x::Point,x₀,R₀,m,ω)
  w = x - x₀
  t = angle(w[1]+w[2]*im)
  w⋅w - (R₀*(1.0+m*sin(ω*t)))^2
end

function popcorn(ranks, parts, nc, ng)
  geo = popcorn()
  pmin = Point(-1.1,-1.1,-1.1)
  pmax = Point(1.1,1.1,1.1)
  bgmodel = CartesianDiscreteModel(ranks,parts,pmin,pmax,(nc,nc,nc);ghost=(ng,ng,ng))
  return bgmodel, geo
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

  t = PArrays.PTimer(ranks)
  PArrays.tic!(t,barrier=true)

  strategy = AggregateCutCellsByThreshold(1.0)
  lcell_to_lroot, lcell_to_root, lcell_to_value =
    map(local_views(cutgeo),cell_indices) do cutgeo,cell_indices
      lid_to_gid = local_to_global(cell_indices)
      aggregate(strategy,cutgeo,geo,lid_to_gid,IN)
    end |> tuple_of_arrays

  PArrays.toc!(t,"New agg - local stage")

  PArrays.tic!(t,barrier=true)

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

  PArrays.toc!(t,"New agg - global stage")

  display(t)

  bgmodel,lcell_to_root
end

function run_benchmark_test(distribute,
                            parts,
                            ncells_x_dir,
                            problem)

  ranks = distribute(LinearIndices((prod(parts),)))
  verbose = i_am_main(ranks)

  verbose && begin
    @info "Parts per direction: $parts"
    @info "Cells in x direction: $ncells_x_dir"
    @info "Ghost layers: $nghost_layers"
    @info "Problem: $(problem==symmetric_kettlebell ? "symmetric" : "asymmetric") kettlebell"
  end
  
  # println("Benchmarking old and new distributed aggregation implementations")
  # println("==============================================")
  # println("Old distributed aggregation")
  # @btime obgmodel,olcell_to_root = run_old_distributed_aggregation(
  #         $ranks,$parts,$ncells_x_dir,$problem)
  # println("New distributed aggregation")
  # @btime nbgmodel,nlcell_to_root = run_new_distributed_aggregation(
  #         $ranks,$parts,$ncells_x_dir,$nghost_layers,$problem)
  # println("==============================================")

  t = PArrays.PTimer(ranks,verbose=true)

  obgmodel,olcell_to_root = run_old_distributed_aggregation(
    ranks,parts,ncells_x_dir,problem)
  for repeat = 1:4
    PArrays.tic!(t,barrier=true)
      obgmodel,olcell_to_root = run_old_distributed_aggregation(
        ranks,parts,ncells_x_dir,problem)
    PArrays.toc!(t,"Old AGG - ncells $ncells_x_dir - run $repeat")
  end

  for nghost_layers in (2,3,4,5)
    nbgmodel,nlcell_to_root = run_new_distributed_aggregation(
      ranks,parts,ncells_x_dir,nghost_layers,problem)
    for repeat = 1:4
      PArrays.tic!(t,barrier=true)
        nbgmodel,nlcell_to_root = run_new_distributed_aggregation(
          ranks,parts,ncells_x_dir,nghost_layers,problem)
      PArrays.toc!(t,"New AGG - ncells $ncells_x_dir - run $repeat - $nghost_layers ghost layers")
    end
  end

  display(t)

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
  #   Triangulation(obgmodel), "data/popcorn_aggregates_old", 
  #   celldata = ["aggregate" => oocell_to_root],
  # );

  # writevtk(
  #   Triangulation(nbgmodel), "data/popcorn_aggregates_new", 
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

function run_distributed_agfem(distribute,
                               parts,
                               ncells_x_dir,
                               nghost_layers,
                               problem)

  ranks = distribute(LinearIndices((prod(parts),)))
  bgmodel, geo = problem(ranks, parts, ncells_x_dir, nghost_layers)
  cutgeo = cut(bgmodel, geo)

  cell_gids = get_cell_gids(bgmodel)
  cell_indices = partition(cell_gids)

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

  ocell_to_root = map(lcell_to_root,own_to_local(cell_gids)) do agg,o_to_l
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
  # TO-DO: Jordi reviews if aggtrian can be created from aggmodel
  # See if swapping the cell_gids of the cutgeo for the new ones
  function active_aggregate_conforming_trian(model,cutgeo,agg_cell_indices)
    trians = map( local_views(cutgeo),
                  agg_cell_indices ) do cutgeo, agg_cell_indices
      Triangulation(with_ghost,agg_cell_indices,cutgeo,ACTIVE)
    end
    aggmodel = GridapDistributed.GenericDistributedDiscreteModel(
      local_views(model),PRange(agg_cell_indices))
    GridapDistributed.DistributedTriangulation(trians,aggmodel)
  end

  aggtrian = active_aggregate_conforming_trian(bgmodel,cutgeo,agg_cell_indices)

  writevtk(aggtrian,"data/kettlebell_aggtrian");

  order = 1
  reffe = ReferenceFE(lagrangian,Float64,order)

  V = TestFESpace(aggtrian,reffe)
  # (TMP) This gives wrong output when root is not owned, but does not affect below.
  lcell_to_lroot = _local_aggregates(lcell_to_root,cell_gids)

  aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs = AgFEMSpace(V,lcell_to_lroot,agg_cell_indices)

  # # Proper sanity check to use when full table of constraints computed
  # Vagg = AgFEMSpace(V,lcell_to_lroot,agg_cell_indices)
  # u(x) = x[1]+x[2]
  # uh = interpolate_everywhere(u,Vagg)
  # dΩ = Measure(aggtrian,2)
  # @test ∑(∫(u-uh)dΩ) < 1.0e-12

  # map(ranks,local_views(spaces)) do r,space
  #   aggtrian = get_triangulation(space)
  #   writevtk(
  #     aggtrian,
  #     "data/kettlebell_aggtrian_$(r)",
  #     celldata=[
  #       "part" => fill(r,num_cells(aggtrian)),
  #       # "dof_ids" => string.(get_cell_dof_ids(space))
  #     ],
  #   );
  # end;

end

end