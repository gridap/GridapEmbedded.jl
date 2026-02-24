using MPI
using PartitionedArrays

include("../eric_dev.jl")
import .DistributedAggregationP4estMeshes as TestModule

if !MPI.Initialized()
  MPI.Init()
end

with_mpi() do distribute 
  TestModule.run(distribute)
end

# MPI.Finalize()