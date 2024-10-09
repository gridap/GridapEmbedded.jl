module SequentialTests

using Test

@time @testset "PoissonSeq" begin include("PoissonTests.jl") end
@time @testset "DiscreteGeoPoissonSeq" begin include("DistributedDiscreteGeometryPoissonTest.jl") end
@time @testset "LSDiscreteGeoPoissonSeq" begin include("DistributedLSDiscreteGeometryPoissonTest.jl") end

## Disabled due to Issue #87
# @time @testset "PeriodicDiscreteGeoPoissonSeq" begin include("PeriodicDistributedDiscreteGeometryPoissonTest.jl") end

end
