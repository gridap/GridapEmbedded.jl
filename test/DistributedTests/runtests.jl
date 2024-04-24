module DistributedTests

if Int === Int64
  include("sequential/runtests.jl")
  include("mpi/runtests.jl")
end

end # module
