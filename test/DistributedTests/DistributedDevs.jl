module DistributedDevs

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

using GridapEmbedded.Distributed: distributed_aggregate
using GridapEmbedded.Distributed: has_remote_aggregation

distribute = PartitionedArrays.DebugArray

np = (2,1)

ranks = distribute(LinearIndices((prod(np),)))

L = 1
p0 = Point(0.0,0.0)
pmin = p0-L/2
pmax = p0+L/2

R = 0.3
x = VectorValue(1.0,0.0)

d1 = VectorValue(0.4,0.0)
d2 = VectorValue(0.0,0.05)
geo = disk(R,x0=p0)
# geo2 = quadrilateral(;x0=p0-d2/2,d1=d1,d2=d2)
# geo = union(geo1,geo2)

n = 9
mesh_partition = (n,n)
bgmodel = CartesianDiscreteModel(ranks,np,pmin,pmax,mesh_partition)

cutgeo = cut(bgmodel,geo)

bgf_to_ioc = compute_bgfacet_to_inoutcut(bgmodel,geo)

Ω = Triangulation(cutgeo)

writevtk(Ω,"trian")

strategy = AggregateCutCellsByThreshold(1.0)
aggregates,aggregate_owner,aggregate_neig = distributed_aggregate(
  strategy,cutgeo,geo,IN)


gids = get_cell_gids(bgmodel)
# TODO: do not use owner for this query (has ghost aggregation?)
local_aggretation = map(aggregate_owner,own_to_local(gids),own_to_owner(gids)) do owns,o_to_l,o
  owners = map(Reindex(owns),o_to_l)
  all(lazy_map(==,owners,o))
end
has_local_aggretation = reduction(&,local_aggretation,destination=:all) |> PartitionedArrays.getany


@test !has_remote_aggregation(bgmodel,aggregates)


@test !has_local_aggretation
# @test !has_remote_aggregation

cut_in = map(aggregates) do agg
  findall(!iszero,agg)
end

cut_owner = map(aggregate_owner,cut_in) do owners,cut_in
  owners[cut_in]
end

aggregates
gids = get_cell_gids(bgmodel)

laggregates = map(aggregates,global_to_local(gids)) do agg,g_to_l
   map(agg) do a
    if iszero(a)
      a
    else
      g_to_l[a]
    end
  end
end

oaggregates = map(aggregates,own_to_local(gids)) do agg,o_to_l
  map(Reindex(agg),o_to_l)
end

oaggregate_owner = map(aggregate_owner,own_to_local(gids)) do agg,o_to_l
  map(Reindex(agg),o_to_l)
end


Ωbg = Triangulation(bgmodel)
Ω = Triangulation(cutgeo)
Ωin = Triangulation(cutgeo,IN)
Γ = EmbeddedBoundary(cutgeo)
Ω_act = Triangulation(cutgeo,ACTIVE)



n_Γ = get_normal_vector(Γ)

order = 1
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)

reffe = ReferenceFE(lagrangian,Float64,order)

Vstd = FESpace(Ω_act,reffe)


map(num_cells,local_views(Ω_act))


map(local_views(Ω_act),ranks) do trian,p
  writevtk(trian,"trian_act_$p")
end

V = AgFEMSpace(bgmodel,Vstd,laggregates)
# V = Vstd
U = TrialFESpace(V)

h = 1.0/n
γd = 10.0

u(x) = x[1] - x[2]
f(x) = -Δ(u)(x)
ud(x) = u(x)


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


writevtk(Ωin,"trian_in")
writevtk(Γ,"bnd")
writevtk(Ωbg,"bgtrian",celldata=
   ["aggregate"=>oaggregates,
    "aggregate_owner"=>oaggregate_owner])

writevtk(Ω,"trian",cellfields=["uh"=>uh])




using GridapEmbedded.Distributed: add_remote_aggregates

np = (3,)

ranks = distribute(LinearIndices((prod(np),)))
m = CartesianDiscreteModel(ranks,np,(0, 1),(9,))


aggregates = DebugArray([[1, 9, 8, 5], [3, 4, 9, 6, 7], [6, 7, 3, 9]])
aggregate_owner = DebugArray([[1, 3, 3, 2], [1, 2, 3, 2, 3], [2, 3, 1, 3]])

am = add_remote_aggregates(m,aggregates,aggregate_owner)

gids = partition(get_cell_gids(am))
test_gids = DebugArray([[1,2,3,4,5,8,9],[3,4,5,6,7,9],[6,7,8,9,3]])
map(gids,test_gids) do gids,_gids
  @test all(gids .==_gids)
end

map(local_views(am),ranks) do m,p
  writevtk(Triangulation(m),"m$p")
end


# import GridapEmbedded.Distributed: change_bgmodel
# using GridapEmbedded.Distributed: DistributedEmbeddedDiscretization
# using GridapDistributed: DistributedDiscreteModel

# TODO:
# - print aggregates [x]
# - move to src [x]
# - test aggregates with several geometries [x]
# - reconstruct paths [x]
# - add remote roots to model [x]
# - migrate model [x]
# - agfem fe space with remotes

# TODO: parallel aggregation (parallel aggfem article)
# 0. exchange measures and aggregates through gid [x]
# 1. root in ghost layer [x]
# 2. root in neighbors [x]
# 3. root in neighbors of neighbors [x]

# TODO:
# - test possion eq with root in ghost layer (simplify this geometry) [x]
# - test graph reconstruction in 1D arrays [x]
# - add remote roots to model [x]

end
