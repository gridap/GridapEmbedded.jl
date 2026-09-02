using MPI
using PartitionedArrays

include("../constraint_merging_distributed.jl")
import .DistributedAgFEM4P4estMeshes as TestModule

if !MPI.Initialized()
  MPI.Init()
end

model, cutgeo, geo, uₕ = with_mpi() do distribute
  # Square:                 begin with initial_uniform_refs = 3
  # SquareWithCircularHole: begin with initial_uniform_refs = 3
  # Flower:                 begin with initial_uniform_refs = 6
  # SinusoidalBand:         begin with initial_uniform_refs = 5
  # el2s, eh1s, hs = run_convergence_test(SquareWithCircularHole(),
  #                      initial_uniform_refs=3,
  #                      order=2,
  #                      num_refinements=4,
  #                      problem=out_fe_space())
  TestModule.run_single_test(distribute,
                      TestModule.Square(),
                      num_parts=4,
                      initial_uniform_refs=3,
                      order=2,
                      problem=TestModule.in_fe_space(order=2),
                      write_solution=true)
end

MPI.Finalize()