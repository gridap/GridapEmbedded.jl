module SequentialTests

using Test

@time @testset "PoissonSeq" begin include("PoissonTests.jl") end
@time @testset "DiscreteGeoPoissonSeq" begin include("DistributedDiscreteGeometryPoissonTest.jl") end
@time @testset "PeriodicPoissonSeq" begin include("PeriodicPoissonTests.jl") end
@time @testset "PeriodicDiscreteGeoPoissonSeq" begin include("PeriodicDistributedDiscreteGeometryPoissonTest.jl") end

end
