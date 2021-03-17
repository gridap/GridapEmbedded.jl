module CellAggregationTests

using Test
using Gridap
using GridapEmbedded
using GridapEmbedded.AgFEM

const R = 0.4
geom = disk(R,x0=Point(0.5,0.5))
n = 30
partition = (n,n)

domain = (0,1,0,1)
bgmodel = CartesianDiscreteModel(domain,partition)

cutdisc = cut(bgmodel,geom)

strategy = AggregateAllCutCells()

aggregates = aggregate(strategy,cutdisc,geom)

colors = color_aggregates(aggregates,bgmodel)

trian = Triangulation(bgmodel)

aggregates_out = aggregate(strategy,cutdisc,geom,OUT)

colors_out = color_aggregates(aggregates_out,bgmodel)

#writevtk(trian,"trian",
#  celldata=["cellin"=>aggregates,"color"=>colors,"cellin_out"=>aggregates_out,"color_out"=>colors_out])

cutdisc_facets = cut_facets(bgmodel,geom)
aggregates = aggregate(strategy,cutdisc,cutdisc_facets)

facet_to_inoutcut = compute_bgfacet_to_inoutcut(bgmodel,geom)

@test aggregates == aggregate(strategy,cutdisc,facet_to_inoutcut)

@test aggregates_out == aggregate(strategy,cutdisc,geom,OUT,facet_to_inoutcut)

end # module
