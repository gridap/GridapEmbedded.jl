module DistributedAggregationMPI

using PartitionedArrays
const PArrays = PartitionedArrays
using MPI

include("../distributed_aggregation.jl")

const DA = DistributedAggregation

if ! MPI.Initialized()
  MPI.Init()
end

problem = DA.symmetric_kettlebell

if MPI.Comm_size(MPI.COMM_WORLD) == 3
  with_mpi() do distribute
    DA.run_benchmark_test(distribute,(3,1),9,2,problem)
    for ncells_x_dir in (9,18,36,72), nghost_layers in (2,3,4)
      DA.run_benchmark_test(distribute,
                         (3,1),
                         ncells_x_dir,
                         nghost_layers,
                         problem)
    end
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 4
  with_mpi() do distribute
    DA.run_benchmark_test(distribute,(2,2),8,2,problem)
    for ncells_x_dir in (8,16,32,64), nghost_layers in (2,3,4)
      DA.run_benchmark_test(distribute,
                         (2,2),
                         ncells_x_dir,
                         nghost_layers,
                         problem)
    end
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 8
  with_mpi() do distribute
    DA.run_benchmark_test(distribute,(4,2),8,2,problem)
    for ncells_x_dir in (16,32,64,128), nghost_layers in (2,3,4,5)
      DA.run_benchmark_test(distribute,
                         (4,2),
                         ncells_x_dir,
                         nghost_layers,
                         problem)
    end
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 12
  with_mpi() do distribute
    DA.run_benchmark_test(distribute,(4,3),12,2,problem)
    for ncells_x_dir in (24,48,96,192), nghost_layers in (2,3,4,5)
      DA.run_benchmark_test(distribute,
                         (4,3),
                         ncells_x_dir,
                         nghost_layers,
                         problem)
    end
  end
end

end # module