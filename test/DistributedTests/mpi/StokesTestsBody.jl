module MPITestsBody

using PartitionedArrays
const PArrays = PartitionedArrays
using MPI

include("../StokesTests.jl")

if ! MPI.Initialized()
  MPI.Init()
end

function all_tests(distribute,parts)
  ranks = distribute(LinearIndices((prod(parts),)))

  t = PArrays.PTimer(ranks,verbose=true)
  PArrays.tic!(t)

  StokesTests.main(distribute,parts,cells=(2,2))
  PArrays.toc!(t,"Stokes")

  display(t)
end

if MPI.Comm_size(MPI.COMM_WORLD) == 16
  with_mpi() do distribute
    all_tests(distribute,(4,4))
  end
elseif MPI.Comm_size(MPI.COMM_WORLD) == 4
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
