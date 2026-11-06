module InterfacesTests

using Test

@testset "SubFaceTriangulations" begin include("SubFaceTriangulationsTests.jl") end

@testset "EmbeddedDiscretizations" begin include("EmbeddedDiscretizationsTests.jl") end

@testset "EmbeddedFacetDiscretizations" begin include("EmbeddedFacetDiscretizationsTests.jl") end

@testset "Cutters" begin include("CuttersTests.jl") end

@testset "ImmersedDiscreteModels" begin include("ImmersedDiscreteModelsTests.jl") end

@testset "Issue Reproducers" begin
  @testset "Issue115" begin include("issue_115.jl") end
end

end # module
