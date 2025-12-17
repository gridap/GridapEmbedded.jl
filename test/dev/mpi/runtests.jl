module DistributedAggregationMPITests

using Test
using MPI

mpidir = @__DIR__
testdir = joinpath(mpidir,"..")
repodir = joinpath(testdir,"..","..")

function run_test(procs,file)
  mpiexec() do cmd
    run(`$cmd -n $procs $(Base.julia_cmd()) --project=$repodir $(joinpath(mpidir,file))`)
  end
end

run_test(12288,"runtests_body.jl")

end # module
