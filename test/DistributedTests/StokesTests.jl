module StokesTests

using Gridap
using Gridap.ReferenceFEs
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

u(x) = VectorValue(x[1]*x[1],x[2])
p(x) = x[1] - x[2]

f(x) = - Δ(u)(x) + ∇(p)(x)
g(x) = (∇⋅u)(x)

function main(distribute,parts;n=4,cells=(n,n),order=2)

  @show parts
  ranks = distribute(LinearIndices((prod(parts),)))

  φ(x) = x[1]*x[1] + x[2]*x[2] - 1.0
  # phi = AlgoimCallLevelSetFunction(
  #   x -> ( x[1]*x[1] + x[2]*x[2] + x[3]*x[3] ) - 1.0,
  #   x -> VectorValue( 2.0 * x[1], 2.0 * x[2], 2.0 * x[3] ) )

  domain = (-1.1,1.1,-1.1,1.1)
  model₀ = CartesianDiscreteModel(ranks,parts,domain,cells)

  orderᵠ = 2
  Ω₀ = Triangulation(model₀)
  reffeᵠ = ReferenceFE(lagrangian,Float64,orderᵠ)
  V = TestFESpace(Ω₀,reffeᵠ,conformity=:H1)
  φₕ = interpolate_everywhere(φ,V)
  phi = AlgoimCallLevelSetFunction(φₕ,∇(φₕ))

  degree = order == 1 ? 3 : 2*order
  vquad = Quadrature(algoim,phi,degree,phase=IN)
  v_cell_quad,is_a = CellQuadratureAndActiveMask(model₀,vquad)
  squad = Quadrature(algoim,phi,degree)
  s_cell_quad,is_c = CellQuadratureAndActiveMask(model₀,squad)
  
  model,aggregates = aggregate(model₀,is_a,is_c,IN)

  Ω = Triangulation(model)
  Ωᵃ,dΩᵃ = TriangulationAndMeasure(Ω,v_cell_quad,is_a)

  Λᵃ = SkeletonTriangulation(model)
  dΛᵃ = Measure(Λᵃ,2*order)

  _,dΓ = TriangulationAndMeasure(Ω,s_cell_quad,is_c)
  n_Γ = normal(phi,Ω)

  # dbg = parts[1]
  # writevtk(dΓ,"cut_quad_$dbg")

  # @show √(∑( ∫( 1.0 )dΩᵃ ))
  # @show √(∑( ∫( 1.0 )dΓ  ))

  reffeᵘ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  reffeˢ = ReferenceFE(lagrangian,VectorValue{2,Float64},order,space=:S)
  reffeᵖ = ReferenceFE(lagrangian,Float64,order-1,space=:P)
  reffeˡ = ReferenceFE(lagrangian,Float64,order-2,space=:P)

  Vstdᵘ = TestFESpace(Ωᵃ,reffeᵘ)
  Vstdᵖ = TestFESpace(Ωᵃ,reffeᵖ)

  Vserᵘ = TestFESpace(Ωᵃ,reffeˢ,conformity=:L2)

  Vᵘ = AgFEMSpace(model,Vstdᵘ,aggregates,Vserᵘ,reffeˢ)
  Vᵖ = AgFEMSpace(model,Vstdᵖ,aggregates)
  Vˡ = Vstdˡ

  Uᵘ = TrialFESpace(Vᵘ)
  Uᵖ = TrialFESpace(Vᵖ)
  Uˡ = TrialFESpace(Vˡ)
  
  Y = MultiFieldFESpace([Vᵘ,Vᵖ,Vˡ])
  X = MultiFieldFESpace([Uᵘ,Uᵖ,Uˡ])

  # Nitsche method
  γᵈ = 2.0*order^2
  τ₁ = 500.0
  h = (domain[2]-domain[1])/cells[1]

  a((u,p,l),(v,q,ℓ)) =
    ∫( ∇(u)⊙∇(v) - q*(∇⋅u) - p*(∇⋅v) )dΩᵃ +
    ∫( p*ℓ )dΩᵃ + ∫( q*l )dΩᵃ +
    ∫( (τ₁/h)*jump(l)*jump(ℓ) )dΛᵃ +
    ∫( (γᵈ/h)*(u⋅v) - 
        v⋅(n_Γ⋅∇(u)) - u⋅(n_Γ⋅∇(v)) + 
        p*(n_Γ⋅v) + q*(n_Γ⋅u) )dΓ

  l((v,q,ℓ)) =
    ∫( v⋅f - q*g )dΩᵃ +
    ∫( (γᵈ/h)*(u⋅v) - u⋅(n_Γ⋅∇(v)) + (q*n_Γ)⋅u )dΓ

  # FE problem
  op = AffineFEOperator(a,l,X,Y)
  uₕ,pₕ,lₕ = solve(op)

  # writevtk(Ωᵃ,"stokes",cellfields=["u"=>uₕ,"p"=>pₕ,"l"=>lₕ])

  euₕ = u - uₕ
  epₕ = p - pₕ

  l2(e) = √(∑( ∫( e⋅e )dΩᵃ ))
  h1(e) = √(∑( ∫( e⋅e + ∇(e)⊙∇(e) )dΩᵃ ))

  el2ᵘ = l2(euₕ)
  eh1ᵘ = h1(euₕ)
  el2ᵖ = l2(epₕ)
  # el2ˡ = l2(lₕ)

  # @test el2ᵘ < 1.0e-014
  # @test eh1ᵘ < 1.0e-013
  # @test el2ᵖ < 1.0e-014
  @show el2ᵘ
  @show eh1ᵘ
  @show el2ᵖ
  # @show el2ˡ

end

function main_zeromean(distribute,parts;n=4,cells=(n,n),order=2)

  @show parts
  ranks = distribute(LinearIndices((prod(parts),)))

  φ(x) = x[1]*x[1] + x[2]*x[2] - 1.0
  # phi = AlgoimCallLevelSetFunction(
  #   x -> ( x[1]*x[1] + x[2]*x[2] + x[3]*x[3] ) - 1.0,
  #   x -> VectorValue( 2.0 * x[1], 2.0 * x[2], 2.0 * x[3] ) )

  domain = (-1.1,1.1,-1.1,1.1)
  model₀ = CartesianDiscreteModel(ranks,parts,domain,cells)

  orderᵠ = 2
  Ω₀ = Triangulation(model₀)
  reffeᵠ = ReferenceFE(lagrangian,Float64,orderᵠ)
  V = TestFESpace(Ω₀,reffeᵠ,conformity=:H1)
  φₕ = interpolate_everywhere(φ,V)
  phi = AlgoimCallLevelSetFunction(φₕ,∇(φₕ))

  degree = order == 1 ? 3 : 2*order
  vquad = Quadrature(algoim,phi,degree,phase=IN)
  v_cell_quad,is_a = CellQuadratureAndActiveMask(model₀,vquad)
  squad = Quadrature(algoim,phi,degree)
  s_cell_quad,is_c = CellQuadratureAndActiveMask(model₀,squad)
  
  model,aggregates = aggregate(model₀,is_a,is_c,IN)
  # map(aggregates) do agg
  #   @show agg
  # end

  Ω = Triangulation(model)
  Ωᵃ,dΩᵃ = TriangulationAndMeasure(Ω,v_cell_quad,is_a)
  _,dΓ = TriangulationAndMeasure(Ω,s_cell_quad,is_c)
  n_Γ = normal(phi,Ω)
  
  reffeᵘ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  reffeᵖ = ReferenceFE(lagrangian,Float64,order-1,space=:P)

  Vstdᵘ = TestFESpace(Ωᵃ,reffeᵘ)
  Vstdᵖ = TestFESpace(Ωᵃ,reffeᵖ,constraint=:zeromean)
  
  Vᵘ = AgFEMSpace(model,Vstdᵘ,aggregates)
  Vᵖ = AgFEMSpace(model,Vstdᵖ,aggregates)

  Uᵘ = TrialFESpace(Vᵘ)
  Uᵖ = TrialFESpace(Vᵖ)
  
  Y = MultiFieldFESpace([Vᵘ,Vᵖ])
  X = MultiFieldFESpace([Uᵘ,Uᵖ])

  # Nitsche method
  γᵈ = 2.0*order^2
  h = (domain[2]-domain[1])/cells[1]

  a((u,p),(v,q)) =
    ∫( ∇(u)⊙∇(v) - q*(∇⋅u) - p*(∇⋅v) )dΩᵃ +
    ∫( (γᵈ/h)*(u⋅v) - 
        v⋅(n_Γ⋅∇(u)) - u⋅(n_Γ⋅∇(v)) + 
        p*(n_Γ⋅v) + q*(n_Γ⋅u) )dΓ

  l((v,q)) =
    ∫( v⋅f - q*g )dΩᵃ +
    ∫( (γᵈ/h)*(u⋅v) - u⋅(n_Γ⋅∇(v)) + (q*n_Γ)⋅u )dΓ

  # FE problem
  op = AffineFEOperator(a,l,X,Y)
  uₕ,pₕ = solve(op)
  # writevtk(Ωᵃ,"stokes",cellfields=["u"=>uₕ,"p"=>pₕ])

  euₕ = u - uₕ
  epₕ = p - pₕ

  l2(e) = √(∑( ∫( e⋅e )dΩᵃ ))
  h1(e) = √(∑( ∫( e⋅e + ∇(e)⊙∇(e) )dΩᵃ ))

  el2ᵘ = l2(euₕ)
  eh1ᵘ = h1(euₕ)
  el2ᵖ = l2(epₕ)

  # @test el2ᵘ < 1.0e-014
  # @test eh1ᵘ < 1.0e-013
  # @test el2ᵖ < 1.0e-014
  @show el2ᵘ
  @show eh1ᵘ
  @show el2ᵖ

end

end # module
