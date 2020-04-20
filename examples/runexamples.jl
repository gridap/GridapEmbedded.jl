module RunExamples

using Test

@time @testset "poisson_csg_cutfem" begin
  include("poisson_csg_cutfem.jl")
  PoissonCSGCutFEM.main(n=30)
end

@time @testset "stokes_tube_with_obstacle_cutfem" begin
  include("stokes_tube_with_obstacle_cutfem.jl")
  StokesTubeWithObstacleCutFEM.main(n=10)
end

end # module
