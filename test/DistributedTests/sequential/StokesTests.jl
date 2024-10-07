module StokesTestsSeq
using PartitionedArrays
include("../StokesTests.jl")
with_debug() do distribute
    StokesTests.main(distribute,(1,1),cells=(16,16))
    StokesTests.main(distribute,(2,2),cells=(16,16))
end
end