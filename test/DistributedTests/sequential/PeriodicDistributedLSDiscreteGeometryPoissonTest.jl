module PeriodicDistributedLSDiscreteGeometryPoissonTestsSeq
using PartitionedArrays
include("../PeriodicDistributedLSDiscreteGeometryPoissonTest.jl")
with_debug() do distribute
  PeriodicDistributedLSDiscreteGeometryPoissonTest.main(distribute,(2,2),cells=(21,21))
end
end
