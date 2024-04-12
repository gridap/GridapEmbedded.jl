module AggregatesTestsSeq
using PartitionedArrays
include("../AggregatesTests.jl")
with_debug() do distribute
    DistributedAggregatesTests.main(distribute,(2,2))
end
end
