using Gridap
using GridapDistributed
using GridapP4est
using PartitionedArrays
using MPI

using GridapEmbedded.AlgoimUtils

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute

  # Initial Octree Distributed Discrete Model
  ranks        = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))
  coarse_model = CartesianDiscreteModel((0,1,0,1),(1,1))
  num_levels_initial_refinement = 2
  dmodel = OctreeDistributedDiscreteModel(ranks,
                                          coarse_model,
                                          num_levels_initial_refinement)

  # Coarsen first four cells
  fmodel_refine_coarsen_flags = 
    map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
      flags = zeros(Int,length(indices))
      flags[1:4] .=  coarsen_flag
      flags
  end
  fmodel,_ = Gridap.Adaptivity.adapt(dmodel,fmodel_refine_coarsen_flags);
  # writevtk(fmodel,"fmodel")

  # Interpolate a fun on a point of the non-conforming grid
  Ω = Triangulation(fmodel)
  order = 1
  reffe_s = ReferenceFE(lagrangian,Float64,order)
  Qₕ = TestFESpace(Ω,reffe_s,conformity=:H1)
  φ = zero(Qₕ)

  sm = Gridap.CellData.KDTreeSearch(num_nearest_vertices=3)
  iφ = map(local_views(φ)) do liφ
    Gridap.CellData.Interpolable(liφ,searchmethod=sm)
  end |> GridapDistributed.DistributedInterpolable
  iφ(Point(0.35,0.35)) # OK

  φ(Point(0.35,0.35)) # KO. Using default num_nearest_vertices = 1 fails because the only candidate vertex is a hanging node and does not have in its list of cells the coarse cell that contains the point.

  true
end

# MPI.Finalize()