module PeriodicPoissonTestsSeq
using PartitionedArrays
include("../PeriodicPoissonTests.jl")
with_debug() do distribute
    PeriodicPoissonTests.main(distribute,(2,2),cells=(21,21))
    PeriodicPoissonTests.main(distribute,(4,1),cells=(21,21),geometry=:remotes)
end
end
