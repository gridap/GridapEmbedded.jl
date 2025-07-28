module LevelSetCuttersRun

using Test

@testset "AnalyticalGeometries" begin include("AnalyticalGeometriesTests.jl") end

@testset "DiscreteGeometries" begin include("DiscreteGeometriesTests.jl") end

@testset "LookupTables" begin include("LookupTablesTests.jl") end

@testset "CutTriangulations" begin include("CutTriangulationsTests.jl") end

@testset "LevelSetCutters" begin include("LevelSetCuttersTests.jl") end

@testset "GeometricalDifferentiation" begin include("GeometricalDifferentiationTests.jl") end

end # module
