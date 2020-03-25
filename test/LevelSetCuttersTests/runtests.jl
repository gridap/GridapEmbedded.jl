module LevelSetCuttersTests

using Test

@testset "LookupTables" begin include("LookupTablesTests.jl") end

@testset "Geometries" begin include("GeometriesTests.jl") end

@testset "SubTriangulations" begin include("SubTriangulationsTests.jl") end

@testset "LevelSetCutters" begin include("LevelSetCuttersTests.jl") end

end # module
