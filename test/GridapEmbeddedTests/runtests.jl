module GridapEmbeddedTests

using Test

@testset "PoissonCutFEM" begin include("PoissonCutFEMTests.jl") end

@testset "BimaterialPoissonCutFEM" begin include("BimaterialPoissonCutFEMTests.jl") end

@testset "EmbeddedBimaterialPoissonCutFEM" begin include("EmbeddedBimaterialPoissonCutFEMTests.jl") end

end
