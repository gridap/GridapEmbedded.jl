module AlgoimTests

using Test

@testset "AlgoimInterface" begin include("AlgoimInterfaceTests.jl") end

@testset "QuadratureDegree" begin include("QuadratureDegreeTests.jl") end

@testset "DualQuadratures" begin include("DualQuadraturesTests.jl") end

@testset "PoissonAlgoim" begin include("PoissonAlgoimTests.jl") end

@testset "ClosestPoint" begin include("ClosestPointTests.jl") end

@testset "VolumeConservation" begin include("VolumeConservationTests.jl") end

@testset "Visualization" begin include("VisualizationTests.jl") end

end # module