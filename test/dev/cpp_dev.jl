
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

  # Coarsen first for cells
  fmodel_refine_coarsen_flags = 
    map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
      flags = zeros(Int,length(indices))
      # flags[1:4] .=  coarsen_flag
      flags
  end
  fmodel,_ = Gridap.Adaptivity.adapt(dmodel,fmodel_refine_coarsen_flags);
  writevtk(fmodel,"fmodel")

  Ω = Triangulation(fmodel)

  order = 3
  reffe_s = ReferenceFE(lagrangian,Float64,order)
  reffe_v = ReferenceFE(lagrangian,VectorValue{2,Float64},order)

  Qₕ = TestFESpace(Ω,reffe_s,conformity=:H1)
  Vₕ = TestFESpace(Ω,reffe_v,conformity=:H1)

  max_refinement_level = num_levels_initial_refinement
  cpps = compute_closest_point_projections(
    fmodel,Vₕ,φ_fun,order,max_refinement_level,cppdegree=3)
  dists = compute_distance_fe_function(
    fmodel,Qₕ,Vₕ,φ_fun,order,max_refinement_level,cppdegree=3)

  writevtk(Ω,"Ω",
    cellfields=["cpp"=>cpps,"dist"=>dists],nsubcells=1)

  # # Interpolable on OctreeDistributedDiscreteModel 
  # # is apparently not implemented yet
  # _φ = interpolate_everywhere(val,Qₕ)
  # φ_field = AlgoimCallLevelSetFunction(_φ,∇(_φ))
  # cpps = compute_closest_point_projections(fmodel,Vₕ,φ_field,max_refinement_level)
  # writevtk(Ω,"Ω",cellfields=["cpp"=>cpps,"phi"=>φ_field.values∘cpps])

  # ISSUES
  # - For order > 2 and nsubcells > 1, wrong cpps and dists, in spite
  #   of the fact that free_dof_values coincide with the sequential case
  # - Can we interpolate anywhere in OctreeDistributedDiscreteModels?

  true
end

# MPI.Finalize()