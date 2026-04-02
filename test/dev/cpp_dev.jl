using Gridap
using GridapEmbedded
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

  # Level set function from an analytical expression
  R = 0.4
  val(x) = (x[1]-0.5)^2 + (x[2]-0.5)^2 - R^2
  grad(x) = ∇(y -> val(y),x)
  φ_fun = AlgoimCallLevelSetFunction(val,grad)

  # # Check against serial
  # model = CartesianDiscreteModel((0,1,0,1),(2,2))
  # order = 3
  # reffe_s = ReferenceFE(lagrangian,Float64,order)
  # Ω = Triangulation(model)
  # V = TestFESpace(Ω,reffe_s,conformity=:H1)
  # dists = compute_distance_fe_function(model,V,φ_fun,order,cppdegree=3)
  # writevtk(Ω,"Ω",cellfields=["dist"=>dists],nsubcells=1)

  # Coarsen first four cells
  fmodel_refine_coarsen_flags = 
    map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
      flags = zeros(Int,length(indices))
      flags[1:4] .=  coarsen_flag
      flags
  end
  fmodel,_ = Gridap.Adaptivity.adapt(dmodel,fmodel_refine_coarsen_flags);
  writevtk(fmodel,"fmodel")

  Ω = Triangulation(fmodel)

  order = 3
  reffe_s = ReferenceFE(lagrangian,Float64,order)
  reffe_v = ReferenceFE(lagrangian,VectorValue{2,Float64},order)

  Qₕ = TestFESpace(Ω,reffe_s,conformity=:H1) # Scalar FE space for distance function
  Vₕ = TestFESpace(Ω,reffe_v,conformity=:H1) # Vector FE space for closest point projections

  max_refinement_level = num_levels_initial_refinement

  # # RMK: Algoim CPP algorithms work on uniform grids
  # # In order to use them on non-uniform grids, we work
  # # on an upper bound grid, corresponding to the uniform 
  # # mesh obtained at maximum refinement level. On this 
  # # maximal grid, the working arrays for the CPP are
  # # computed using the coordinates of the grid points 
  # # and its level set values.
  # cpps = compute_closest_point_projections(
  #   fmodel,Vₕ,φ_fun,order,max_refinement_level,cppdegree=3)
  # dists = compute_distance_fe_function(
  #   fmodel,Qₕ,Vₕ,φ_fun,order,max_refinement_level,cppdegree=3)

  # writevtk(Ω,"Ω",cellfields=["cpp"=>cpps,"dist"=>dists],nsubcells=3)

  _φ = interpolate_everywhere(val,Qₕ)
  φ_field = AlgoimCallLevelSetFunction(_φ,∇(_φ))
  sm = Gridap.CellData.KDTreeSearch(num_nearest_vertices=3)
  iφ_field = map(local_views(φ_field.values)) do iφ
    Gridap.CellData.Interpolable(iφ,searchmethod=sm)
  end |> GridapDistributed.DistributedInterpolable
  cpps = compute_closest_point_projections(
    fmodel,Vₕ,φ_field,order,max_refinement_level)
  writevtk(Ω,"Ω",cellfields=["cpp"=>cpps,"phi"=>iφ_field∘cpps],nsubcells=3)

  # [TODO]: Compute distance function
  # [TODO]: Check visualisation of cpp and distance function

  true
end

# MPI.Finalize()