module MPITestsBody

using PartitionedArrays
const PArrays = PartitionedArrays
using MPI

include("../PoissonTests.jl")
include("../AggregatesTests.jl")

if ! MPI.Initialized()
  MPI.Init()
end

function all_tests(distribute,parts)
  ranks = distribute(LinearIndices((prod(parts),)))

  t = PArrays.PTimer(ranks,verbose=true)
  PArrays.tic!(t)

  PoissonTests.main(distribute,parts)
  PoissonTests.main(distribute,(prod(parts),1),cells=(12,12),geometry=:remotes)
  PArrays.toc!(t,"Poisson")

  if prod(parts) == 4
    DistributedAggregatesTests.main(distribute,parts)
  end
  PArrays.toc!(t,"Aggregates")

  display(t)
end

if MPI.Comm_size(MPI.COMM_WORLD) == 4
  with_mpi() do distribute
    all_tests(distribute,(2,2))
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 1
  with_mpi() do distribute
    all_tests(distribute,(1,1))
  end
else
  MPI.Abort(MPI.COMM_WORLD,0)
end

end #module
