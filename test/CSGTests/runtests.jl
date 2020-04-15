module CSGTests

using Test

@testset "Nodes" begin include("NodesTests.jl") end

@testset "Geometries" begin include("GeometriesTests.jl") end

end # module
