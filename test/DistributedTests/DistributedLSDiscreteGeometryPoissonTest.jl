module DistributedLSDiscreteGeometryPoissonTest

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

using GridapEmbedded.CSG
using GridapEmbedded.LevelSetCutters

function main(distribute,parts;
  threshold=1,
  n=21,
  cells=(n,n))

  order = 2
  ranks = distribute(LinearIndices((prod(parts),)))
  domain = (0,1,0,1)
  bgmodel = CartesianDiscreteModel(ranks,parts,domain,cells)

  u(x) = (x[1]-0.5)^2 + (x[2]-0.5)^2
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  reffe = ReferenceFE(lagrangian,Float64,order)
  Ω_bg = Triangulation(bgmodel)
  V_bg = FESpace(Ω_bg,reffe)
  φh = interpolate(x->sqrt((x[1]-0.5)^2+(x[2]-0.5)^2)-0.35,V_bg)
  geo = DiscreteGeometry(φh,bgmodel)

  D = 2
  cell_meas = map(get_cell_measure∘Triangulation,local_views(bgmodel))
  meas = map(first,cell_meas) |> PartitionedArrays.getany
  h = meas^(1/D)

  cutgeo = cut(bgmodel,geo)
  cutgeo_facets = cut_facets(bgmodel,geo)

  strategy = AggregateCutCellsByThreshold(threshold)
  bgmodel,cutgeo,aggregates = aggregate(strategy,cutgeo)

  Ω_act = Triangulation(cutgeo,ACTIVE)
  Ω = Triangulation(cutgeo,PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)

  n_Γ = get_normal_vector(Γ)

  order = 2
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓ = Measure(Γ,degree)

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

  path = mktempdir()
  writevtk(Ω_bg,joinpath(path,"trian"),
    celldata=[
      "aggregate"=>own_aggregates,
      "color"=>own_colors,
      "gid"=>own_to_global(gids)])#,
  #  cellfields=["uh"=>uh])

  writevtk(Ω,joinpath(path,"trian_O"),cellfields=["uh"=>uh])
  writevtk(Γ,joinpath(path,"trian_G"))
  @test el2/ul2 < 1.e-8
  @test eh1/uh1 < 1.e-7

end

end # module