module GridapEmbeddedTests

using Test

@time @testset "PoissonCutFEM" begin include("PoissonCutFEMTests.jl") end

@time @testset "PoissonAgFEM" begin include("PoissonAgFEMTests.jl") end

@time @testset "BimaterialPoissonCutFEM" begin include("BimaterialPoissonCutFEMTests.jl") end

@time @testset "EmbeddedBimaterialPoissonCutFEM" begin include("EmbeddedBimaterialPoissonCutFEMTests.jl") end

@time @testset "StokesCutFEM" begin include("StokesCutFEMTests.jl") end

@time @testset "StokesAgFEM" begin include("StokesAgFEMTests.jl") end

@time @testset "TraceFEM" begin include("TraceFEMTests.jl") end

end
