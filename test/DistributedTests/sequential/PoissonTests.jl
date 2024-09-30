module PoissonTestsSeq
using PartitionedArrays
include("../PoissonTests.jl")
with_debug() do distribute
    # PoissonTests.main(distribute,(2,2))
    # PoissonTests.main(distribute,(1,4),cells=(4,8))
    # PoissonTests.main(distribute,(2,4),cells=(8,8))
    # PoissonTests.main(distribute,(4,1),cells=(12,12),geometry=:remotes)
    # PoissonTests.main_algoim(distribute,(2,1,2),cells=(4,4,4),solution_degree=2)
    # PoissonTests.main_algoim(distribute,(4,1,2),cells=(8,8,8),solution_degree=2)
    # PoissonTests.main_algoim(distribute,(1,1,4),cells=(8,8,8),solution_degree=2)
    # PoissonTests.main_algoim(distribute,(2,2,2),cells=(8,8,8))
    PoissonTests.main_algoim(distribute,(2,2),cells=(4,4))
    PoissonTests.main_algoim(distribute,(4,4),cells=(16,16))
end
end