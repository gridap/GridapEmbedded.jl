module PeriodicDistributedDiscreteGeometryPoissonTest
using PartitionedArrays
include("../PeriodicDistributedDiscreteGeometryPoissonTest.jl")
with_debug() do distribute
  PeriodicDistributedDiscreteGeometryPoissonTest.main(distribute,(4,4),cells=(31,31),geometry=:circle)
  PeriodicDistributedDiscreteGeometryPoissonTest.main(distribute,(4,1),cells=(31,31),geometry=:remotes)
end
end
