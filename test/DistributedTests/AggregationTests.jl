module DistributedAggregationTests

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

using GridapEmbedded.Distributed: distributed_aggregate

distribute = PartitionedArrays.DebugArray

np = (4,1)

ranks = distribute(LinearIndices((prod(np),)))

L = 1
p0 = Point(0.0,0.0)
pmin = p0-L/2
pmax = p0+L/2

R = 0.2
x = VectorValue(1.0,0.0)

d1 = VectorValue(0.8,0.0)
d2 = VectorValue(0.0,0.02)
geo1 = disk(R,x0=p0-x*(L/2-R))
geo2 = quadrilateral(;x0=p0-d1/2-d2/2,d1=d1,d2=d2)
geo = union(geo1,geo2)

n = 16
mesh_partition = (n,n)
bgmodel = CartesianDiscreteModel(ranks,np,pmin,pmax,mesh_partition;ghost=(2,2))

cutgeo = cut(bgmodel,geo)

bgf_to_ioc = compute_bgfacet_to_inoutcut(bgmodel,geo)

Ω = Triangulation(cutgeo)

strategy = AggregateCutCellsByThreshold(1.0)
aggregates,aggregate_owner,aggregate_neig = distributed_aggregate(
  strategy,cutgeo,geo,IN)

cut_in = map(aggregates) do agg
  findall(!iszero,agg)
end

cut_owner = map(aggregate_owner,cut_in) do owners,cut_in
  owners[cut_in]
end

map(ranks,cut_owner) do p,cut_owner
  if p ∈ (3,4)
    @test all(i->i==2,cut_owner)
  end
end

gids = get_cell_gids(bgmodel)

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

path = mktempdir()
writevtk(Ωin,joinpath(path,"trian_in"))
writevtk(Γ,joinpath(path,"bnd"))
writevtk(Ωbg,joinpath(path,"bgtrian"),celldata=
   ["aggregate"=>oaggregates,
    "aggregate_owner"=>oaggregate_owner])


# TODO:
# - print aggregates [x]
# - move to src [x]
# - test aggregates with several geometries [x]
# - reconstruct paths
# - add remote roots to model

# TODO: parallel aggregation (parallel aggfem article)
# 0. exchange measures and aggregates through gid
# 1. root in ghost layer
# 2. root in neighbors
# 3. root in neighbors of neighbors


# TODO:
# - test possion eq with root in ghost layer (simplify this geometry)
# - test graph reconstruction in 1D arrays
# - add remote roots to model

end # module
