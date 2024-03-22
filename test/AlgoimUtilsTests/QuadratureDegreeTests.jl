module QuadratureDegreeTests

using Test
using Gridap
using Gridap.ReferenceFEs
using GridapEmbedded
using GridapEmbedded.Interfaces

φ(x) = x[1]-x[2]
∇φ(x) = VectorValue(1.0,-1.0,0.0)
phi = AlgoimCallLevelSetFunction(φ,∇φ)

function run_diffusion_reaction_bulk(domain,cells,order,degree)

  u(x) = sum(x)^order
  f(x) = u(x)-Δ(u)(x)
  ud(x) = u(x)

  model = CartesianDiscreteModel(domain,cells)
  Ω = Triangulation(model)

  vquad = Quadrature(algoim,phi,degree,phase=IN)
  Ωᵃ,dΩᵃ,is_a = TriangulationAndMeasure(Ω,vquad)
  squad = Quadrature(algoim,phi,degree)
  _,dΓ,is_c = TriangulationAndMeasure(Ω,squad)
  n_Γ = normal(phi,Ω)

  cell_to_cellin = aggregate(Ω,is_a,is_c,IN)

  reffe = ReferenceFE(lagrangian,Float64,order)
  Vstd = TestFESpace(Ωᵃ,reffe,conformity=:H1,dirichlet_tags="boundary")

  V = AgFEMSpace(Vstd,cell_to_cellin)
  U = TrialFESpace(V,ud)

  # Nitsche method
  γᵈ = 5.0*order^2
  h = (domain[2]-domain[1])/cells[1]

  a(u,v) =
    ∫( v⋅u )dΩᵃ +
    ∫( ∇(v)⋅∇(u) )dΩᵃ +
    ∫( (γᵈ/h)*v*u - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )dΓ

  l(v) =
    ∫( v*f )dΩᵃ +
    ∫( (γᵈ/h)*v*ud - (n_Γ⋅∇(v))*ud )dΓ

  # FE problem
  op = AffineFEOperator(a,l,U,V)
  uₕ = solve(op)

  eₕ = u - uₕ

  l2(u) = √(∑( ∫( u*u )dΩᵃ ))
  h1(u) = √(∑( ∫( u*u + ∇(u)⋅∇(u) )dΩᵃ ))

  el2 = l2(eₕ)
  eh1 = h1(eₕ)

  el2, eh1

end

left_project(u,n) = u - n⊗(n⋅u)
∇ᵈ(u,n) = left_project(∇(u),n)

function run_diffusion_reaction_surface(domain,cells,order,degree)

  u(x) = sum(x)^order
  f(x) = u(x)-Δ(u)(x)
  ud(x) = u(x)

  model = CartesianDiscreteModel(domain,cells)
  Ω = Triangulation(model)

  squad = Quadrature(algoim,phi,degree)
  Ωc,dΓ = TriangulationAndMeasure(Ω,squad)

  n_Γ = normal(phi,Ωc)
  dΩc = Measure(Ωc,2*order)

  reffe = ReferenceFE(lagrangian,Float64,order)

  V = TestFESpace(Ωc,reffe,dirichlet_tags="boundary")
  U = TrialFESpace(V,ud)

  m(u,v) = ∫( v⋅u )dΓ
  a(u,v) = ∫( ∇ᵈ(v,n_Γ)⋅∇ᵈ(u,n_Γ) )dΓ

  # Γ-orthogonal derivative volume stabilisation
  h = (domain[2]-domain[1])/cells[1]
  γˢ = 10.0*h
  s(u,v) = ∫( γˢ*((n_Γ⋅∇(v))⊙(n_Γ⋅∇(u))) )dΩc

  aₛ(u,v) = m(u,v) + a(u,v) + s(u,v)
  l(v) = ∫( v⊙f )dΓ

  op = AffineFEOperator(aₛ,l,U,V)
  uₕ = solve(op)
  eₕ = uₕ - u

  l2(v) = √(∑(∫(v⋅v)dΓ))
  l2(eₕ),l2(∇ᵈ(eₕ,n_Γ))

end

domain = (-1.1,1.1,-1.1,1.1,-1.1,1.1)
partition = (8,8,8)

# BULK tests

order = 1
degree = 2
el2, eh1 = run_diffusion_reaction_bulk(domain,partition,order,degree)
@test el2 < 1.0e-015
@test eh1 < 1.0e-014

order = 2
degree = 3
el2, eh1 = run_diffusion_reaction_bulk(domain,partition,order,degree)
@test el2 < 1.0e-013
@test eh1 < 1.0e-012

# For order > 2, expect accuracy loss with Lagrangian AgFEM

# order = 3
# degree = 5
# el2, eh1 = run_diffusion_reaction_bulk(domain,partition,order,degree)
# @test el2 < 1.0e-011
# @test eh1 < 1.0e-009

# order = 4
# degree = 6
# el2, eh1 = run_diffusion_reaction_bulk(domain,partition,order,degree)
# @test el2 < 1.0e-007
# @test eh1 < 1.0e-005

# SURFACE tests

order = 1
degree = 2
el2, eh1 = run_diffusion_reaction_surface(domain,partition,order,degree)
@test el2 < 1.0e-015
@test eh1 < 1.0e-014

order = 2
degree = 3
el2, eh1 = run_diffusion_reaction_surface(domain,partition,order,degree)
@test el2 < 1.0e-013
@test eh1 < 1.0e-012

# For order > 2, expect accuracy loss with Lagrangian AgFEM

# order = 3
# degree = 5
# el2, eh1 = run_diffusion_reaction_surface(domain,partition,order,degree)
# @test el2 < 1.0e-009
# @test eh1 < 1.0e-008

# order = 4
# degree = 6
# el2, eh1 = run_diffusion_reaction_surface(domain,partition,order,degree)
# @test el2 < 1.0e-007
# @test eh1 < 1.0e-005

end # module