module CellAggregationSpaceTests

using Gridap
using GridapEmbedded
using GridapEmbedded.LevelSetCutters
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Test

# Background Mesh
domain = (0,2,0,1)
n_x = 10
n_t = 10
partition = (n_x,n_t)
bgmodel = CartesianDiscreteModel(domain,partition)

# domain
geo = quadrilateral(x0=Point(0,0),d1=VectorValue(1,0),d2=VectorValue(1,1))
cutgeo = cut(bgmodel,geo)
model = DiscreteModel(cutgeo)

#AgFEM Strategy
strategy = AggregateSpaceCutCells()
aggregates = aggregatespace(strategy,cutgeo)

colors = color_aggregates_space(aggregates,bgmodel)

trian = Triangulation(bgmodel)

#d = "./"
#fi = joinpath(d,"trian")
#writevtk(trian,fi,celldata=["colors"=>colors,"aggregates"=>aggregates])

end  # module CellAggregationSpaceTests
