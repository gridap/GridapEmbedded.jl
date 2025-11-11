module DistributedAggregationP4estMeshes

  using Gridap
  using GridapEmbedded
  using GridapEmbedded.Interfaces: CUT
  using GridapDistributed
  using PartitionedArrays
  using MPI

  using Gridap.Geometry: get_faces

  using P4est_wrapper
  using GridapP4est

  include("NonConformingGridTopologies.jl")
  using .NonConformingGridTopologies

  function run(distribute)

    ranks = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))
    coarse_model = CartesianDiscreteModel((-1,1,-1,1),(1,1))
    num_uniform_refinements = 1
    num_ghost_layers = 1
    D = num_dims(coarse_model)

    dmodel = OctreeDistributedDiscreteModel(ranks,
                                            coarse_model,
                                            num_uniform_refinements;
                                            num_ghost_layers=num_ghost_layers)
    
    geo1 = quadrilateral(;x0=Point(-0.9,-0.9),
                         d1=VectorValue(1.8,0.0),
                         d2=VectorValue(0.0,1.8))
    geo2 = ! disk(0.4)
    geo = intersect(geo1,geo2)

    cutgeo = cut(dmodel, geo)
    cell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geo)
    fmodel_refine_coarsen_flags = 
      map(ranks,
          partition(get_cell_gids(dmodel.dmodel)),
          cell_to_inoutcut) do rank,indices,cell_to_inoutcut
        flags = zeros(Int,length(indices))
        flags .= nothing_flag
        toref = findall(c->c==CUT,cell_to_inoutcut)
        flags[toref] .= refine_flag
        flags
    end
    fmodel,_ = Gridap.Adaptivity.adapt(dmodel,fmodel_refine_coarsen_flags);

    for i in 1:5
      cutgeo = cut(fmodel, geo)
      cell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geo)
      fmodel_refine_coarsen_flags = 
        map(ranks,
            partition(get_cell_gids(fmodel)),
            cell_to_inoutcut) do rank,indices,cell_to_inoutcut
          flags = zeros(Int,length(indices))
          flags .= nothing_flag
          toref = findall(c->c==CUT,cell_to_inoutcut)
          flags[toref] .= refine_flag
          flags
      end
      fmodel,_ = Gridap.Adaptivity.adapt(fmodel,fmodel_refine_coarsen_flags);
    end
    fmodel, _ = GridapDistributed.redistribute(fmodel)

    cutgeo = cut(fmodel, geo)
    writevtk(EmbeddedBoundary(cutgeo),"data/quad_bnd");
    # The next line fails unless I comment out the assertion 
    # in the constructor of AppendedTriangulations.jl:
    # @assert get_background_model(a) === get_background_model(b)
    writevtk(Triangulation(cutgeo,ACTIVE_IN),"data/quad_phys");
    # writevtk(Triangulation(cutgeo,PHYSICAL_IN),"data/quad_phys");

    cell_gids = get_cell_gids(fmodel)
    cell_indices = partition(cell_gids)

    ncgt = NonConformingGridTopology(fmodel)
    strategy = AggregateCutCellsByThreshold(1.0)
    lcell_to_lroot, lcell_to_root, lcell_to_value =
      map(local_views(cutgeo),cell_indices,ncgt) do cutgeo,cell_indices,ncgt
        lid_to_gid = local_to_global(cell_indices)
        aggregate(strategy,cutgeo,geo,lid_to_gid,IN,grid_topology=ncgt)
      end |> tuple_of_arrays

    # Output for verification of lcell_to_root map
    ocell_to_root = map(lcell_to_root,own_to_local(cell_gids)) do agg,o_to_l
      map(Reindex(agg),o_to_l)
    end

    writevtk(EmbeddedBoundary(cutgeo),"data/quad_bnd");
    writevtk(Triangulation(cutgeo,PHYSICAL),"data/quad_phys");
    writevtk(
      Triangulation(fmodel), "data/quad_aggregates", 
      celldata = ["aggregate" => ocell_to_root],
    );

    # map(ncgt) do ncgt
    #   get_faces(ncgt,D,0)
    #   get_faces(ncgt,0,D)
    #   get_faces(ncgt,D,1)
    #   get_faces(ncgt,1,D)
    #   get_faces(ncgt,D,2)
    #   get_faces(ncgt,2,D)
    #   # get_faces(ncgt,D,3)
    #   # get_faces(ncgt,3,D)
    # end

    return nothing

  end

end