module StokesTestsSeq
using PartitionedArrays
include("../StokesTests.jl")
with_debug() do distribute
    StokesTests.main(distribute,(1,1),cells=(2,2),name="a")
    # StokesTests.main(distribute,(1,2),cells=(2,2),name="b")
    # StokesTests.main(distribute,(2,1),cells=(2,2),name="c")
    StokesTests.main(distribute,(2,2),cells=(2,2),name="d")
    # StokesTests.main_zeromean(distribute,(1,1),cells=(4,4))
    # StokesTests.main_zeromean(distribute,(1,2),cells=(4,4))
    # StokesTests.main_zeromean(distribute,(2,1),cells=(4,4))
    # StokesTests.main_zeromean(distribute,(2,2),cells=(4,4))
    # StokesTests.main(distribute,(4,4),cells=(16,16))
end
end