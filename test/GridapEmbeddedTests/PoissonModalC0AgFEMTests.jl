module PoissonModalC0AgFEMTests

  using Gridap
  using GridapEmbedded
  using Test

  using GridapEmbedded.Interfaces: CUT, IN

  u(x) = (x[1]+x[2])^3
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  const R = 0.42

  n = 6
  order = 3

  geom = disk(R,x0=Point(0.5,0.5))
  partition = (n,n)
  domain = (0,1,0,1)

  bgmodel = CartesianDiscreteModel(domain,partition)
  h = (domain[2]-domain[1])/n
  cutgeo = cut(bgmodel,geom)

  Ω_act = Triangulation(cutgeo,ACTIVE)
  model = get_active_model(Ω_act)

  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo,geom)
  bboxes = compute_cell_to_dface_bboxes(model,cutgeo,aggregates)

  Ω_bg = Triangulation(bgmodel)
  Ω = Triangulation(cutgeo,PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)

  n_Γ = get_normal_vector(Γ)

  cutdeg, degree = 2*num_dims(model)*order, 2*order
  dΩ = Measure(Ω,cutdeg,degree)
  dΓ = Measure(Γ,cutdeg)

  # reffe = ReferenceFE(lagrangian,Float64,order)
  reffe = ReferenceFE(modalC0,Float64,order,bboxes)
  Vstd = TestFESpace(Ω_act,reffe,conformity=:H1)
  V = AgFEMSpace(Vstd,aggregates) # lagrangian
  # V = AgFEMSpace(Vstd,aggregates,modalC0)
  U = TrialFESpace(V)

  γd = 5.0*order^2

  a(u,v) =
    ∫( ∇(v)⋅∇(u) )dΩ +
    ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )dΓ

  l(v) =
    ∫( v*f )dΩ +
    ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud )dΓ

  op = AffineFEOperator(a,l,U,V)

  uh = solve(op)
 
  e = u - uh

  l2(u) = sqrt(abs(sum( ∫( u*u )dΩ )))
  h1(u) = sqrt(abs(sum( ∫( u*u + ∇(u)⋅∇(u) )dΩ )))

  el2 = l2(e)
  eh1 = h1(e)
  ul2 = l2(uh)
  uh1 = h1(uh)

  @test el2/ul2 < 1.e-8
  @test eh1/uh1 < 1.e-7

end # module