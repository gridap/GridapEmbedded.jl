module AlgoimTests

using Test

@time @testset "AlgoimInterface" begin include("AlgoimInterfaceTests.jl") end

@time @testset "QuadratureDegree" begin include("QuadratureDegreeTests.jl") end

@time @testset "DualQuadratures" begin include("DualQuadraturesTests.jl") end

@time @testset "PoissonAlgoim" begin include("PoissonAlgoimTests.jl") end

@time @testset "ClosestPoint" begin include("ClosestPointTests.jl") end

@time @testset "VolumeConservation" begin include("VolumeConservationTests.jl") end

@time @testset "Visualization" begin include("VisualizationTests.jl") end

end # module