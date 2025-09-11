module DistributedAggregationMPITests

using Test
using MPI

mpidir = @__DIR__
testdir = mpidir
repodir = joinpath(testdir,"..","..")

function run_test(procs,file)
  mpiexec() do cmd
    run(`$cmd -n $procs $(Base.julia_cmd()) --project=$repodir $(joinpath(mpidir,file))`)
  end
end

run_test( 3,"distributed_aggregation.jl")
run_test( 4,"distributed_aggregation.jl")
run_test( 8,"distributed_aggregation.jl")
run_test(12,"distributed_aggregation.jl")

end # module
