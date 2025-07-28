module DistributedLSDiscreteGeometryPoissonTestsSeq
using PartitionedArrays
include("../DistributedLSDiscreteGeometryPoissonTest.jl")
with_debug() do distribute
  DistributedLSDiscreteGeometryPoissonTest.main(distribute,(2,2),cells=(21,21))
end
end
