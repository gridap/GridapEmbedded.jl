module SequentialTests

using Test

@time @testset "PoissonSeq" begin include("PoissonTests.jl") end
@time @testset "PoissonSeq" begin include("DistributedDiscreteGeometryPoissonTest.jl") end

end
