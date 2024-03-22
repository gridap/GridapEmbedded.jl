module AlgoimInterfaceTests

using Test
using Gridap
using Gridap.ReferenceFEs
using GridapEmbedded
using GridapEmbedded.Interfaces

function run_case_raw_point(degree,phase)

  u(x) = (x[1]*x[1]/(1.0*1.0)+x[2]*x[2]/(0.5*0.5)) - 1.0
  gradu(x) = VectorValue(2.0*x[1]/(1.0*1.0),2.0*x[2]/(0.5*0.5))

  phi = AlgoimCallLevelSetFunction(u,gradu)

  xmin = Point(-1.1,-1.1)
  xmax = Point(1.1,1.1)

  _, weights = fill_quad_data(phi,xmin,xmax,phase,degree)

  s = ∑(weights)
  if phase == IN
    @test s ≈ 1.6055602335693198
  elseif phase == CUT
    @test s ≈ 5.04463087731472
  else
    error()
  end

end

function run_case_function(cells,degree,phase)

  domain = (-1.1,1.1,-1.1,1.1)

  model = CartesianDiscreteModel(domain,cells)
  Ω = Triangulation(model)

  u(x) = (x[1]*x[1]/(1.0*1.0)+x[2]*x[2]/(0.5*0.5)) - 1.0
  gradu(x) = VectorValue(2.0*x[1]/(1.0*1.0),2.0*x[2]/(0.5*0.5))

  phi = AlgoimCallLevelSetFunction(u,gradu)
  quad = Quadrature(algoim,phi,degree,phase=phase)
  dΩ = Measure(Ω,quad,data_domain_style=PhysicalDomain())

  s = ∫(1)dΩ
  if phase == IN
    @test ∑(s) ≈ 1.570796326817732
  elseif phase == CUT
    @test ∑(s) ≈ 4.844224109684756
  else
    error()
  end

end

function run_case_fe_function(cells,order,degree,phase)

  domain = (-1.1,1.1,-1.1,1.1,-1.1,1.1)

  model = CartesianDiscreteModel(domain,cells)
  Ω = Triangulation(model)

  reffe = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω,reffe)
  U = TrialFESpace(V)

  u(x) = (x[1]*x[1]/(1.0*1.0)+x[2]*x[2]/(0.5*0.5)+x[3]*x[3]/(0.33*0.33)) - 1.0
  gradu(x) = VectorValue(2.0*x[1]/(1.0*1.0),2.0*x[2]/(0.5*0.5),2.0*x[3]/(0.33*0.33))

  uₕ = interpolate(u,U)
  ∇uₕ = ∇(uₕ)

  phi = AlgoimCallLevelSetFunction(uₕ,∇uₕ)
  quad = Quadrature(algoim,phi,degree,phase=phase)
  dΩ = Measure(Ω,quad,data_domain_style=PhysicalDomain())

  s = ∫(1)dΩ
  if phase == IN
    @test ∑(s) ≈ 0.6911505144863921
  elseif phase == CUT
    @test ∑(s) ≈ 4.382658365204242
  else
    error()
  end

end

f(x) = (x[1]*x[1]/(1.0*1.0)+x[2]*x[2]/(0.5*0.5)+x[3]*x[3]/(0.33*0.33)) - 1.0
gradf(x) = VectorValue(2.0*x[1]/(1.0*1.0),2.0*x[2]/(0.5*0.5),2.0*x[3]/(0.33*0.33))

phi = AlgoimCallLevelSetFunction(f,gradf)
@test phi(VectorValue(1.0,0.5,0.33)) == 2.0
@test all(phi.∇φ(VectorValue(1.0,0.5,0.33)) .≈ VectorValue(2.0, 4.0, 6.0606060606060606))

run_case_raw_point(3,IN)
run_case_raw_point(3,CUT)

run_case_function((64,64),3,IN)
run_case_function((64,66),3,IN)

run_case_function((64,64),3,CUT)
run_case_function((64,66),3,CUT)

run_case_fe_function((32,32,32),2,3,IN)
run_case_fe_function((32,32,32),2,3,CUT)

end # module