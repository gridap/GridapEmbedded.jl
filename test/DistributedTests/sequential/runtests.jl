module SequentialTests

using Test

@time @testset "PoissonSeq" begin include("PoissonTests.jl") end

end
