module MomentFittingTests

using Test

@testset "CutCellMoments" begin include("CutCellMomentsTests.jl") end

@testset "JacobiPolynomialBases" begin include("JacobiPolynomialBasesTests.jl") end

end # module
