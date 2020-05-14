module AgFEMTests

using Test

@testset "CellAggregation" begin include("CellAggregationTests.jl") end

@testset "AgFEMSpaces" begin include("AgFEMSpacesTests.jl") end

end # module
