module DistributedDevs

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

using GridapEmbedded.CSG

## Parallel aggregation

using GridapEmbedded.Distributed: DistributedEmbeddedDiscretization

function _aggregate(
  stragegy,
  cutgeo::DistributedEmbeddedDiscretization,
  in_or_out::Integer=IN)

  map(local_views(cutgeo)) do cutgeo
    aggregate(stragegy,cutgeo,cutgeo.geo,in_or_out)
  end
end

function _aggregate(
  stragegy,
  cutgeo::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out::Integer=IN)

  bgf_to_ioc = compute_bgfacet_to_inoutcut(bgmodel,geo)
  _aggregate(stragegy,cutgeo,geo,in_or_out,bgf_to_ioc)
end

function _aggregate(
  stragegy,
  cutgeo::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out::Integer,
  bgf_to_ioc::AbstractArray{<:AbstractVector})

  map(local_views(cutgeo),bgf_to_ioc) do cutgeo,bgf_to_ioc
    #input cell meas or IN cells
    aggregate(stragegy,cutgeo,geo,in_or_out,bgf_to_ioc)
  end
end

# Driver

distribute = identity

np = (2,2)

ranks = distribute(LinearIndices((prod(np),)))

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

n = 8
mesh_partition = (n,n)
bgmodel = CartesianDiscreteModel(ranks,np,pmin,pmax,mesh_partition)
bgtrian = Triangulation(bgmodel)
writevtk(bgtrian,"bgtrian")

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

using GridapEmbedded.Interfaces: CUT
Ω_in = Triangulation(cutgeo,IN)
Ω_out = Triangulation(cutgeo,OUT)
Ω_cut = Triangulation(cutgeo,CUT)

writevtk(Ω_in,"trian_in")
writevtk(Ω_out,"trian_out")
writevtk(Ω_cut,"trian_cut")



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
# V = Vstd
U = TrialFESpace(V)


const γd = 10.0

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

##

distribute = identity

np = (2,2)

ranks = distribute(LinearIndices((prod(np),)))

L = 1
p0 = Point(0.0,0.0)
pmin = p0-L/2
pmax = p0+L/2


R = 0.35
geo = disk(R,x0=p0)

n = 8
mesh_partition = (n,n)
bgmodel = CartesianDiscreteModel(ranks,np,pmin,pmax,mesh_partition)

cutgeo = cut(bgmodel,geo)

bgf_to_ioc = compute_bgfacet_to_inoutcut(bgmodel,geo)

# BUG: aggregates do not have information of ghost cut cells
# TODO: exchange measures and aggregates
# _aggregate(strategy,cutgeo)


# TODO: parallel aggregation (parallel aggfem article)
# 0. exchange measures and aggregates through gid
# 1. root in ghost layer
# 2. root in neighbors
# 3. root in neighbors of neighbors

end # module
