module DistributedAggregationMPI

using PartitionedArrays
const PArrays = PartitionedArrays
using MPI

include("../distributed_aggregation.jl")

const DA = DistributedAggregation

if ! MPI.Initialized()
  MPI.Init()
end

problem = DA.disk

if MPI.Comm_size(MPI.COMM_WORLD) == 4
  with_mpi() do distribute
    DA.run_benchmark_test(distribute,(2,2),8,2,problem)
    for ncells_x_dir in (8,16,32,64)
      DA.run_benchmark_test(distribute,
                         (2,2),
                         ncells_x_dir,
                         problem)
    end
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 24
  with_mpi() do distribute
    DA.run_benchmark_test(distribute,(4,6),24,2,problem)
    for ncells_x_dir in (96,192,384,768,1536)
      DA.run_benchmark_test(distribute,
                         (4,6),
                         ncells_x_dir,
                         problem)
    end
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 96
  with_mpi() do distribute
    DA.run_benchmark_test(distribute,(8,12),48,2,problem)
    for ncells_x_dir in (192,384,768,1536,3072)
      DA.run_benchmark_test(distribute,
                         (8,12),
                         ncells_x_dir,
                         problem)
    end
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 384
  with_mpi() do distribute
    DA.run_benchmark_test(distribute,(16,24),96,2,problem)
    for ncells_x_dir in (384,768,1536,3072,6144)
      DA.run_benchmark_test(distribute,
                         (16,24),
                         ncells_x_dir,
                         problem)
    end
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 1536
  with_mpi() do distribute
    DA.run_benchmark_test(distribute,(32,48),192,2,problem)
    for ncells_x_dir in (768,1536,3072,6144,12288)
      DA.run_benchmark_test(distribute,
                         (32,48),
                         ncells_x_dir,
                         problem)
    end
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 6144
  with_mpi() do distribute
    DA.run_benchmark_test(distribute,(64,96),384,2,problem)
    for ncells_x_dir in (1536,3072,6144,12288)
      DA.run_benchmark_test(distribute,
                         (64,96),
                         ncells_x_dir,
                         problem)
    end
  end
end

# problem = DA.popcorn

# if MPI.Comm_size(MPI.COMM_WORLD) == 24
#   with_mpi() do distribute
#     DA.run_benchmark_test(distribute,(4,3,2),24,2,problem)
#     for ncells_x_dir in (24,48,96)
#       DA.run_benchmark_test(distribute,
#                          (4,3,2),
#                          ncells_x_dir,
#                          problem)
#     end
#   end
# elseif MPI.Comm_size(MPI.COMM_WORLD) == 192
#   with_mpi() do distribute
#     DA.run_benchmark_test(distribute,(8,6,4),192,2,problem)
#     for ncells_x_dir in (48,96,192)
#       DA.run_benchmark_test(distribute,
#                          (8,6,4),
#                          ncells_x_dir,
#                          problem)
#     end
#   end
# elseif MPI.Comm_size(MPI.COMM_WORLD) == 1536
#   with_mpi() do distribute
#     DA.run_benchmark_test(distribute,(16,12,8),96,2,problem)
#     for ncells_x_dir in (96,192,384)
#       DA.run_benchmark_test(distribute,
#                          (16,12,8),
#                          ncells_x_dir,
#                          problem)
#     end
#   end
# elseif MPI.Comm_size(MPI.COMM_WORLD) == 12288
#   with_mpi() do distribute
#     DA.run_benchmark_test(distribute,(32,24,16),192,2,problem)
#     for ncells_x_dir in (192,384,768)
#       DA.run_benchmark_test(distribute,
#                          (32,24,16),
#                          ncells_x_dir,
#                          problem)
#     end
#   end
# end

end # module