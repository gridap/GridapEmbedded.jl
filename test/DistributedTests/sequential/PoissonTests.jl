module PoissonTestsSeq
using PartitionedArrays
include("../PoissonTests.jl")
with_debug() do distribute
    PoissonTests.main(distribute,(2,2))
    PoissonTests.main(distribute,(4,1),cells=(12,12),geometry=:remotes)
end
end
