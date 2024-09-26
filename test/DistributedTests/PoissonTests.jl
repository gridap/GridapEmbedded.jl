module PoissonTests

using Gridap
using Gridap.ReferenceFEs
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

using GridapEmbedded.CSG

function main(distribute,parts;
  threshold=1,
  n=8,
  cells=(n,n),
  geometry=:circle)

  ranks = distribute(LinearIndices((prod(parts),)))

  u(x) = x[1] - x[2]
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  geometries = Dict(
    :circle => circle_geometry,
    :remotes => remotes_geometry,
  )

  bgmodel,geo = geometries[geometry](ranks,parts,cells)

  D = 2
  cell_meas = map(get_cell_measure∘Triangulation,local_views(bgmodel))
  meas = map(first,cell_meas) |> PartitionedArrays.getany
  h = meas^(1/D)

  cutgeo = cut(bgmodel,geo)
  cutgeo_facets = cut_facets(bgmodel,geo)

  strategy = AggregateCutCellsByThreshold(threshold)
  bgmodel,cutgeo,aggregates = aggregate(strategy,cutgeo)

  Ω_bg = Triangulation(bgmodel)
  Ω_act = Triangulation(cutgeo,ACTIVE)
  Ω = Triangulation(cutgeo,PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)

  n_Γ = get_normal_vector(Γ)

  order = 1
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓ = Measure(Γ,degree)

  reffe = ReferenceFE(lagrangian,Float64,order)

  Vstd = FESpace(Ω_act,reffe)

  V = AgFEMSpace(bgmodel,Vstd,aggregates)
  U = TrialFESpace(V)


  γd = 10.0

  a(u,v) =
    ∫( ∇(v)⋅∇(u) ) * dΩ +
    ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ

  l(v) =
    ∫( v*f ) * dΩ +
    ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

  op = AffineFEOperator(a,l,U,V)
  uh = solve(op)

  e = u - uh

  l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
  h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

  el2 = l2(e)
  eh1 = h1(e)
  ul2 = l2(uh)
  uh1 = h1(uh)

  #
  colors = map(color_aggregates,aggregates,local_views(bgmodel))
  gids = get_cell_gids(bgmodel)


  global_aggregates = map(aggregates,local_to_global(gids)) do agg,gid
    map(i-> i==0 ? 0 : gid[i],agg)
  end
  own_aggregates = map(global_aggregates,own_to_local(gids)) do agg,oid
    map(Reindex(agg),oid)
  end
  own_colors = map(colors,own_to_local(gids)) do col,oid
    map(Reindex(col),oid)
  end

  writevtk(Ω_bg,"trian",
    celldata=[
      "aggregate"=>own_aggregates,
      "color"=>own_colors,
      "gid"=>own_to_global(gids)])#,
  #  cellfields=["uh"=>uh])

  writevtk(Ω,"trian_O",cellfields=["uh"=>uh,"eh"=>e])
  writevtk(Γ,"trian_G")
  @test el2/ul2 < 1.e-8
  @test eh1/uh1 < 1.e-7

end

function main_algoim(distribute,parts;
  threshold=1,
  n=4,
  cells=(n,n,n),
  order=1,
  solution_degree=1)

  ranks = distribute(LinearIndices((prod(parts),)))

  φ(x) = x[1]-x[2]+0.01
  ∇φ(x) = VectorValue(1.0,-1.0,0.0)
  # phi = AlgoimCallLevelSetFunction(
  #   x -> ( x[1]*x[1] + x[2]*x[2] + x[3]*x[3] ) - 1.0,
  #   x -> VectorValue( 2.0 * x[1], 2.0 * x[2], 2.0 * x[3] ) )
 
  u(x) = sum(x)^solution_degree
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  domain = (-1.1,1.1,-1.1,1.1,-1.1,1.1)
  model = CartesianDiscreteModel(ranks,parts,domain,cells)

  Ω₀ = Triangulation(model)
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω₀,reffe,conformity=:H1)
  φₕ = interpolate_everywhere(φ,V)
  phi = AlgoimCallLevelSetFunction(φₕ,∇(φₕ))

  degree = order == 1 ? 3 : 2*order
  vquad = Quadrature(algoim,phi,degree,phase=IN)
  v_cell_quad,is_a = CellQuadratureAndActiveMask(model,vquad)
  squad = Quadrature(algoim,phi,degree)
  s_cell_quad,is_c = CellQuadratureAndActiveMask(model,squad)
  
  model,aggregates = aggregate(model,is_a,is_c,IN)
  # map(aggregates) do agg
  #   @show agg
  # end

  Ω = Triangulation(model)
  Ωᵃ,dΩᵃ = TriangulationAndMeasure(Ω,v_cell_quad,is_a)
  _,dΓ = TriangulationAndMeasure(Ω,s_cell_quad,is_c)

  n_Γ = normal(phi,Ω)

  # writevtk(model,"bgmodel")
  # writevtk(Ω,"trian")
  # writevtk(Ωᵃ,"atrian")

  # colors = map(color_aggregates,aggregates,local_views(model))
  # gids = get_cell_gids(model)

  # global_aggregates = map(aggregates,local_to_global(gids)) do agg,gid
  #   map(i-> i==0 ? 0 : gid[i],agg)
  # end
  # own_aggregates = map(global_aggregates,own_to_local(gids)) do agg,oid
  #   map(Reindex(agg),oid)
  # end
  # own_colors = map(colors,own_to_local(gids)) do col,oid
  #   map(Reindex(col),oid)
  # end

  # writevtk(Ω,"dtrian",
  #   celldata=[
  #     "aggregate"=>own_aggregates,
  #     "color"=>own_colors,
  #     "gid"=>own_to_global(gids)])
  
  Vstd = TestFESpace(Ωᵃ,reffe,conformity=:H1,dirichlet_tags="boundary")
  
  V = AgFEMSpace(model,Vstd,aggregates)
  U = TrialFESpace(V,ud)

  # eΩᵃ = GridapDistributed.add_ghost_cells(Ωᵃ)
  # cell_dof_ids = map(get_cell_dof_ids,local_views(V),local_views(eΩᵃ))
  # map(cell_dof_ids) do cell_dof_ids
  #   @show cell_dof_ids
  # end

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

  # writevtk(Ωᵃ,"rtrian",cellfields=["uh"=>uₕ,"eh"=>eₕ])

  l2(u) = √(∑( ∫( u*u )dΩᵃ ))
  h1(u) = √(∑( ∫( u*u + ∇(u)⋅∇(u) )dΩᵃ ))

  el2 = l2(eₕ)
  eh1 = h1(eₕ)

  # @test el2 < 1.0e-014
  # @test eh1 < 1.0e-013
  @show el2
  @show eh1

end

function circle_geometry(ranks,parts,cells)
  L = 1
  p0 = Point(0.0,0.0)
  pmin = p0-L/2
  pmax = p0+L/2
  R = 0.35
  geo = disk(R,x0=p0)
  bgmodel = CartesianDiscreteModel(ranks,parts,pmin,pmax,cells)
  bgmodel,geo
end

function remotes_geometry(ranks,parts,cells)
  x0 = Point(0.05,0.05)
  d1 = VectorValue(0.9,0.0)
  d2 = VectorValue(0.0,0.1)
  geo1 = quadrilateral(;x0=x0,d1=d1,d2=d2)

  x0 = Point(0.15,0.1)
  d1 = VectorValue(0.25,0.0)
  d2 = VectorValue(0.0,0.6)
  geo2 = quadrilateral(;x0=x0,d1=d1,d2=d2)
  geo = union(geo1,geo2)

  domain = (0, 1, 0, 1)
  bgmodel = CartesianDiscreteModel(ranks,parts,domain,cells)
  bgmodel,geo
end




end # module
