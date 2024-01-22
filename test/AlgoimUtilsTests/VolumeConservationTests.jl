module VolumeConservationTests

using Test
using CxxWrap
using Gridap
using GridapEmbedded

using Gridap.ReferenceFEs

const IN = -1
const OUT = 1
const CUT = 0

domain = (-1.1,1.1,-1.1,1.1)

f(x) = (x[1]*x[1]/(0.5*0.5)+x[2]*x[2]/(0.5*0.5)) - 1.0
function gradf(x::V) where {V}
    V([2.0*x[1]/(0.5*0.5),2.0*x[2]/(0.5*0.5)])
end

function run_case(n::Int,order::Int,cppdegree::Int)

  partition = Int32[n,n]
  bgmodel = CartesianDiscreteModel(domain,partition)
  Ω = Triangulation(bgmodel)

  reffeᵠ = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω,reffeᵠ)

  fₕ = interpolate_everywhere(f,V)
  phi = AlgoimCallLevelSetFunction(fₕ,∇(fₕ))

  degree = order == 1 ? 3 : 2*order
  vquad = Quadrature(algoim,phi,degree,phase=IN)
  dΩ = Measure(Ω,vquad,data_domain_style=PhysicalDomain())
  vol_1 = ∑(∫(1)dΩ)
  # @show vol_1

  cps = compute_closest_point_projections(V,phi,order,cppdegree=cppdegree)

  reffeᶜ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  W = TestFESpace(Ω,reffeᶜ)

  g(x) = VectorValue(1.0,0.0)
  gₕ = interpolate_everywhere(g,W)

  dt = 0.2
  _phi₂ = compute_normal_displacement(cps,phi,gₕ,dt,Ω)
  _phi₂ = get_free_dof_values(fₕ) - _phi₂
  _phi₂ = FEFunction(V,_phi₂)
  phi₂ = AlgoimCallLevelSetFunction(_phi₂,∇(_phi₂))

  vquad = Quadrature(algoim,phi₂,degree,phase=IN)
  dΩ₂ = Measure(Ω,vquad,data_domain_style=PhysicalDomain())
  # writevtk(dΩ₂,"res_vol")
  vol_2 = ∑(∫(1)dΩ₂)
  # @show vol_2

  abs(vol_1-vol_2)/abs(vol_1)

end

order = 2 # Affects CP accuracy _and_ volume conservation
cppdegree = 2 # Affects CP accuracy, but not volume conservation

@test run_case(12,order,cppdegree) < 1.0e-5
# @test run_case(24,order,cppdegree) < 1.0e-6

end # module