module DistributedDiscreteGeometryPoissonTestsSeq
using PartitionedArrays
include("../DistributedDiscreteGeometryPoissonTest.jl")
with_debug() do distribute
  DistributedDiscreteGeometryPoissonTest.main(distribute,(2,2))
  DistributedDiscreteGeometryPoissonTest.main(distribute,(4,1),cells=(12,12),geometry=:remotes)
end
end
