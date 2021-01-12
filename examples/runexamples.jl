module RunExamples

using Test

@time @testset "PoissonCSGCutFEM" begin
  include("PoissonCSGCutFEM/PoissonCSGCutFEM.jl")
  PoissonCSGCutFEM.main(n=30)
end

@time @testset "StokesTubeWithObstacleCutFEM" begin
  include("StokesTubeWithObstacleCutFEM/StokesTubeWithObstacleCutFEM.jl")
  StokesTubeWithObstacleCutFEM.main(n=10)
end

@time @testset "BimaterialLinElastCutFEM" begin
  include("BimaterialLinElastCutFEM/BimaterialLinElastCutFEM.jl")
  BimaterialLinElastCutFEM.main(n=4,outputfile="bimaterial")
end

end # module
