module RunTests

using Test

using MiniQhull

if hasproperty(MiniQhull,:QHULL_WRAPPER_LOADED)
  QHULL_LOADED = MiniQhull.QHULL_WRAPPER_LOADED[]
else
  QHULL_LOADED = true
end

if QHULL_LOADED
  @time @testset "CSG" begin include("CSGTests/runtests.jl") end

  @time @testset "Interfaces" begin include("InterfacesTests/runtests.jl") end

  @time @testset "LevelSetCutters" begin include("LevelSetCuttersTests/runtests.jl") end

  @time @testset "AgFEM" begin include("AgFEMTests/runtests.jl") end

  @time @testset "MomentFitting" begin include("MomentFittingTests/runtests.jl") end

  include(joinpath(@__DIR__,"..","examples","runexamples.jl"))

else

  @warn "MiniQhull not properly installed. Tests are not executed."

end

end # module
