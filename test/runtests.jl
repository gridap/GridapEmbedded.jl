module RunTests

using Test

using MiniQhull

if MiniQhull.QHULL_WRAPPER_LOADED[]

  @time @testset "CSG" begin include("CSGTests/runtests.jl") end
  
  @time @testset "Interfaces" begin include("InterfacesTests/runtests.jl") end
  
  @time @testset "LevelSetCutters" begin include("LevelSetCuttersTests/runtests.jl") end
  
  @time @testset "GridapEmbedded" begin include("GridapEmbeddedTests/runtests.jl") end
  
  include(joinpath(@__DIR__,"..","examples","runexamples.jl"))

else
  @warn "MiniQhull not properly installed. Tests are not executed."

end

end # module
