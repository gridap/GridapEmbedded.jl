module DistributedAggregationP4estMeshes

  using Gridap
  using GridapEmbedded
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
    num_ghost_layers = 2
    D = num_dims(coarse_model)

    dmodel = OctreeDistributedDiscreteModel(ranks,
                                            coarse_model,
                                            num_uniform_refinements;
                                            num_ghost_layers=num_ghost_layers)
    
    fmodel_refine_coarsen_flags = 
      map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
        flags = zeros(Int,length(indices))
        flags .= nothing_flag        
        flags[1] = refine_flag
        flags
    end
    fmodel,_ = Gridap.Adaptivity.adapt(dmodel,fmodel_refine_coarsen_flags);

    ncgt = NonConformingGridTopology(fmodel)
    map(ncgt) do ncgt
      get_faces(ncgt,D,0)
      get_faces(ncgt,0,D)
      get_faces(ncgt,D,1)
      get_faces(ncgt,1,D)
      get_faces(ncgt,D,2)
      get_faces(ncgt,2,D)
      # get_faces(ncgt,D,3)
      # get_faces(ncgt,3,D)
    end

  end

end