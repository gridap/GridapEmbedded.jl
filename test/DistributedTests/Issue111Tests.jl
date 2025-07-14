module DistributedAggregationTests

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

using GridapEmbedded.Distributed: distributed_aggregate

distribute = PartitionedArrays.DebugArray

np = (2,1)

ranks = distribute(LinearIndices((prod(np),)))

n = 8
mesh_partition = (n,n)

L = 1
p0 = Point(0.0,0.0)
pmin = Point(-L/4,-L/8)
pmax = Point(L/4,3*L/8)

R = 0.2
x = VectorValue(1.0,0.0)
h = L/(2*n)

geo1 = disk(R,x0=p0)
geo2 = quadrilateral(;x0=Point(-L/4+h/4,5*h/4),
                      d1=VectorValue(6*h+h/2,0.0),
                      d2=VectorValue(0.0,4*h+h/2))
geo3 = quadrilateral(;x0=Point(-L/4+3*h/4,7*h/4),
                      d1=VectorValue(6*h-h/2,0.0),
                      d2=VectorValue(0.0,4*h-h/2))
geo = union(geo1,intersect(geo2,!geo3))

bgmodel = CartesianDiscreteModel(ranks,np,pmin,pmax,mesh_partition;ghost=(2,2))

cutgeo = cut(bgmodel,geo)

bgf_to_ioc = compute_bgfacet_to_inoutcut(bgmodel,geo)

Ω = Triangulation(cutgeo)

strategy = AggregateCutCellsByThreshold(1.0)
aggregates,aggregate_owner,aggregate_neig = distributed_aggregate(
  strategy,cutgeo,geo,IN)

gids = get_cell_gids(bgmodel)

oaggregates = map(aggregates,own_to_local(gids)) do agg,o_to_l
  map(Reindex(agg),o_to_l)
end

oaggregate_owner = map(aggregate_owner,own_to_local(gids)) do agg,o_to_l
  map(Reindex(agg),o_to_l)
end

map(ranks,oaggregates) do p,a
  if p ∈ (1)
    @test a[end] == 30
  end
end

Ωbg = Triangulation(bgmodel)
Ω = Triangulation(cutgeo)
Ωin = Triangulation(cutgeo,IN)
Γ = EmbeddedBoundary(cutgeo)

path = pwd()
writevtk(Ωin,joinpath(path,"trian_in"))
writevtk(Γ,joinpath(path,"bnd"))
writevtk(Ωbg,joinpath(path,"bgtrian"),celldata=
   ["aggregate"=>oaggregates,
    "aggregate_owner"=>oaggregate_owner])

end # module
