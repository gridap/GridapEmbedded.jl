module DistributedAggregationSeqTests
using PartitionedArrays
include("../distributed_aggregation.jl")
const DA = DistributedAggregation
problem = DA.symmetric_kettlebell
with_debug() do distribute
  DA.run_distributed_agfem(distribute,(3,1),9,2,problem)
end
end # module
