module RunTests

using Test

using Delaunay

  @time @testset "CSG" begin include("CSGTests/runtests.jl") end

  @time @testset "Interfaces" begin include("InterfacesTests/runtests.jl") end

  @time @testset "LevelSetCutters" begin include("LevelSetCuttersTests/runtests.jl") end

  @time @testset "AgFEM" begin include("AgFEMTests/runtests.jl") end

  include(joinpath(@__DIR__,"..","examples","runexamples.jl"))


end # module
