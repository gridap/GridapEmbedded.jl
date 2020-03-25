module GridapEmbeddedTests

using GridapEmbedded
using Test

@time @testset "EmbeddedDiscretizations" begin "EmbeddedDiscretizations.jl" end

@time @testset "Cutters" begin "CuttersTests.jl" end

@time @testset "LookupTables" begin "LookupTablesTests.jl" end

@time @testset "SubTriangulations" begin "SubTriangulationsTests.jl" end

@time @testset "LevelSetCutters" begin "LevelSetCuttersTests.jl" end

end # module
