module PoissonAlgoimTests

using Test
using Gridap
using Gridap.ReferenceFEs
using GridapEmbedded
using GridapEmbedded.Interfaces

function run_poisson(domain,cells,order,solution_degree)

  φ(x) = x[1]-x[2]-1.0
  ∇φ(x) = VectorValue(1.0,-1.0)
  phi = AlgoimCallLevelSetFunction(φ,∇φ)

  u(x) = sum(x)^solution_degree
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  model = CartesianDiscreteModel(domain,cells)
  Ω = Triangulation(model)

  degree = order == 1 ? 3 : 2*order
  vquad = Quadrature(algoim,phi,degree,phase=IN)
  Ωᵃ,dΩᵃ,is_a = TriangulationAndMeasure(Ω,vquad)
  squad = Quadrature(algoim,phi,degree)
  _,dΓ,is_c = TriangulationAndMeasure(Ω,squad)
  n_Γ = normal(phi,Ω)

  aggregates = aggregate(model,is_a,is_c,IN)
  @show aggregates

  reffe = ReferenceFE(lagrangian,Float64,order)
  Vstd = TestFESpace(Ωᵃ,reffe,conformity=:H1,dirichlet_tags="boundary")

  V = AgFEMSpace(Vstd,aggregates)
  U = TrialFESpace(V,ud)

  cell_to_dofs = get_cell_dof_ids(Vstd)
  @show cell_to_dofs

  cell_to_dofs = get_cell_dof_ids(V)
  @show cell_to_dofs

  # Nitsche method
  γᵈ = 2.0*order^2
  h = (domain[2]-domain[1])/cells[1]

  a(u,v) =
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

  @test el2 < 1.0e-014
  @test eh1 < 1.0e-013
  
end

run_poisson((-1.1,1.1,-1.1,1.1),(6,6),1,1)

end # module