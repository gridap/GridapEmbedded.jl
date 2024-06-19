module PeriodicDistributedDiscreteGeometryPoissonTestsSeq
using PartitionedArrays
include("../PeriodicDistributedDiscreteGeometryPoissonTest.jl")
with_debug() do distribute
  PeriodicDistributedDiscreteGeometryPoissonTest.main(distribute,(2,2),cells=(21,21))
  PeriodicDistributedDiscreteGeometryPoissonTest.main(distribute,(4,1),cells=(21,21),geometry=:remotes)
end
end
