module InterfacesTests

using Test

@testset "EmbeddedDiscretizations" begin include("EmbeddedDiscretizationsTests.jl") end

@testset "Cutters" begin include("CuttersTests.jl") end

end # module
