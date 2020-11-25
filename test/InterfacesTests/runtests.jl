module InterfacesTests

using Test

@testset "SubCellTriangulations" begin include("SubCellTriangulationsTests.jl") end

@testset "SubFacetTriangulations" begin include("SubFacetTriangulationsTests.jl") end

@testset "EmbeddedDiscretizations" begin include("EmbeddedDiscretizationsTests.jl") end

@testset "EmbeddedFacetDiscretizations" begin include("EmbeddedFacetDiscretizationsTests.jl") end

@testset "Cutters" begin include("CuttersTests.jl") end

end # module
