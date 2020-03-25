module GridapEmbeddedTests

using Test

@time @testset "Interfaces" begin include("InterfacesTests/runtests.jl") end

@time @testset "LevelSetCutters" begin include("LevelSetCuttersTests/runtests.jl") end

end # module
