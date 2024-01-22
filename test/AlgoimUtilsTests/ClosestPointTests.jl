module ClosestPointTests

using Test
using CxxWrap
using Gridap
using GridapEmbedded

const IN = -1
const OUT = 1
const CUT = 0

function run_case(degree)

  xmin = Point(-1.1,-1.1,-1.1)
  xmax = Point(1.1,1.1,1.1)
  partition = Int32[16,16,16]
  
  domain = (-1.1,1.1,-1.1,1.1,-1.1,1.1)
  bgmodel = CartesianDiscreteModel(domain,partition)
  Ω = Triangulation(bgmodel)
  
  reffeᵠ = ReferenceFE(lagrangian,Float64,1)
  V = TestFESpace(Ω,reffeᵠ,conformity=:L2)
  
  f(x) = (x[1]*x[1]/(0.5*0.5)+x[2]*x[2]/(0.5*0.5)+x[3]*x[3]/(0.5*0.5)) - 1.0
  function gradf(x::V) where {V}
    V([2.0*x[1]/(0.5*0.5),2.0*x[2]/(0.5*0.5),2.0*x[3]/(0.5*0.5)])
  end
  
  fₕ = interpolate_everywhere(f,V)
  phi = AlgoimCallLevelSetFunction(fₕ,∇(fₕ))
  
  coords = fill_cpp_data(phi,partition,xmin,xmax,degree)
  @test maximum(f.(coords)) < 1.0e-4

end

run_case(2)
run_case(3)
run_case(4)
run_case(5)
run_case(-1)

end # module