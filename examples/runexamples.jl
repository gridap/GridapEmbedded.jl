module RunExamples

using Test

@time @testset "stokes_tube_with_obstacle_cutfem" begin
  include("stokes_tube_with_obstacle_cutfem.jl")
  StokesTubeWithObstacleCutFEM.main(n=10)
end

end # module
