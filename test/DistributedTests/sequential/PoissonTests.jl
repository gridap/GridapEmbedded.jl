module PoissonTestsSeq
using PartitionedArrays
include("../PoissonTests.jl")
with_debug() do distribute
    PoissonTests.main(distribute,(2,2))
    PoissonTests.main(distribute,(1,4),cells=(4,8))
    PoissonTests.main(distribute,(2,4),cells=(8,8))
    PoissonTests.main(distribute,(4,1),cells=(12,12),geometry=:remotes)
end
end
