module PeriodicPoissonTests

using Gridap
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

  u(x) = (x[1]-0.5)^2 + (x[2]-0.5)^2
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

  order = 2
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

function circle_geometry(ranks,parts,cells)
  L = 1
  p0 = Point(0.5,0.5)
  pmin = p0-L/2
  pmax = p0+L/2
  R = 0.55
  geo = disk(R,x0=p0)
  bgmodel = CartesianDiscreteModel(ranks,parts,pmin,pmax,cells;isperiodic=(true,true))
  bgmodel,geo
end

function remotes_geometry(ranks,parts,cells)
  x0 = Point(0.0,0.4)
  d1 = VectorValue(1.0,0.0)
  d2 = VectorValue(0.0,0.2)
  geo1 = quadrilateral(;x0=x0,d1=d1,d2=d2)

  x0 = Point(0.4,0.0)
  d1 = VectorValue(0.2,0.0)
  d2 = VectorValue(0.0,1.0)
  geo2 = quadrilateral(;x0=x0,d1=d1,d2=d2)
  geo = union(geo1,geo2)

  domain = (0, 1, 0, 1)
  bgmodel = CartesianDiscreteModel(ranks,parts,domain,cells;isperiodic=(true,true))
  bgmodel,geo
end

end # module
