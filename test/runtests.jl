module RunTests

using Test

@time @testset "Interfaces" begin include("InterfacesTests/runtests.jl") end

@time @testset "LevelSetCutters" begin include("LevelSetCuttersTests/runtests.jl") end

@time @testset "GridapEmbedded" begin include("GridapEmbeddedTests/runtests.jl") end

end # module
