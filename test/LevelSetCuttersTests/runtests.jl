module LevelSetCuttersRun

using Test

@testset "AnalyticalGeometries" begin include("AnalyticalGeometriesTests.jl") end

@testset "DiscreteGeometries" begin include("DiscreteGeometriesTests.jl") end

@testset "LookupTables" begin include("LookupTablesTests.jl") end

@testset "SubTriangulations" begin include("SubTriangulationsTests.jl") end

@testset "LevelSetCutters" begin include("LevelSetCuttersTests.jl") end

end # module
