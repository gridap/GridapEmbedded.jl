module PoissonTests

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

using GridapEmbedded.CSG

function main(distribute,parts;
  n=8,
  cells=(n,n))

  ranks = distribute(LinearIndices((prod(parts),)))

  u(x) = x[1] - x[2]
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  L = 1
  p0 = Point(0.0,0.0)
  pmin = p0-L/2
  pmax = p0+L/2


  R = 0.35
  geo = disk(R,x0=p0)

  R = 0.15
  d = L/4
  geo1 = disk(R,x0=p0-d)
  geo2 = disk(R,x0=p0+d)
  #geo = !union(geo1,geo2)

  bgmodel = CartesianDiscreteModel(ranks,parts,pmin,pmax,cells)
  # bgtrian = Triangulation(bgmodel)
  # writevtk(bgtrian,"bgtrian")

  dp = pmax - pmin
  h = dp[1]/n

  cutgeo = cut(bgmodel,geo)

  strategy = AggregateAllCutCells()
  strategy = AggregateCutCellsByThreshold(0.5)
  aggregates = aggregate(strategy,cutgeo)


  Ω_bg = Triangulation(bgmodel)
  Ω_act = Triangulation(cutgeo,ACTIVE)
  Ω = Triangulation(cutgeo,PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)

  # using GridapEmbedded.Interfaces: CUT
  # Ω_in = Triangulation(cutgeo,IN)
  # Ω_out = Triangulation(cutgeo,OUT)
  # Ω_cut = Triangulation(cutgeo,CUT)

  # writevtk(Ω_in,"trian_in")
  # writevtk(Ω_out,"trian_out")
  # writevtk(Ω_cut,"trian_cut")


  writevtk(Ω_bg,"trian")
  writevtk(Ω_act,"trian_act")
  writevtk(Ω,"trian_O")
  writevtk(Γ,"trian_G")

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

  writevtk(Ω,"trian_O",cellfields=["uh"=>uh])
  writevtk(Γ,"trian_G")
  @test el2/ul2 < 1.e-8
  @test eh1/uh1 < 1.e-7

end

end # module
