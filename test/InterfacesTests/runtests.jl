module InterfacesTests

using Test

@testset "SubTriangulations" begin include("SubTriangulationsTests.jl") end

@testset "FacetSubTriangulations" begin include("FacetSubTriangulationsTests.jl") end

@testset "EmbeddedDiscretizations" begin include("EmbeddedDiscretizationsTests.jl") end

@testset "Cutters" begin include("CuttersTests.jl") end

end # module
