module GridapEmbeddedTests

using Test

@testset "PoissonCutFEM" begin include("PoissonCutFEMTests.jl") end

#@testset "BimaterialPoissonCutFEM" begin include("BimaterialPoissonCutFEMTests.jl") end

end
