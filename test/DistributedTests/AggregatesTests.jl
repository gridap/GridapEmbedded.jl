module DistributedAggregatesTests

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test
using GridapEmbedded.Distributed: distributed_aggregate
using GridapEmbedded.Distributed: has_remote_aggregation
using GridapEmbedded.Distributed: add_remote_aggregates
using GridapEmbedded.Distributed: change_bgmodel
using GridapEmbedded.Distributed: _extract_remote_cells


function main(distribute,parts;vtk=false)

  @assert prod(parts) == 4
  ranks = distribute(LinearIndices((4,)))

  np = (4,1)
  nc = (12,12)

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
  bgmodel = CartesianDiscreteModel(ranks,np,domain,nc)
  cutgeo = cut(bgmodel,geo)
  strategy = AggregateCutCellsByThreshold(1.0)
  aggregates,aggregate_owner = distributed_aggregate(strategy,cutgeo)

  gids = get_cell_gids(bgmodel)
  remote_cells,remote_parts = _extract_remote_cells(gids,aggregates,aggregate_owner)
  _bgmodel = add_remote_aggregates(bgmodel,aggregates,aggregate_owner)
  _cutgeo = change_bgmodel(cutgeo,bgmodel)
  _aggregates = change_bgmodel(aggregates,get_cell_gids(bgmodel))

  @test has_remote_aggregation(bgmodel,aggregates)
  @test ! has_remote_aggregation(_bgmodel,_aggregates)

  _remote_cells = distribute([[],[],[28],[28]])
  map(remote_cells,_remote_cells) do a,b
      @test a == b
  end
  _remote_parts = distribute([[],[],[2],[2]])
  map(remote_parts,_remote_parts) do a,b
      @test a == b
  end

  if vtk
    gids = get_cell_gids(bgmodel)
    Ωbg = Triangulation(bgmodel)
    Γ = EmbeddedBoundary(cutgeo)
    Ω = Triangulation(cutgeo)
    own_aggregates = map(aggregates,own_to_local(gids)) do agg,o_to_l
      agg[o_to_l]
    end
    writevtk(Ωbg,"bgmodel",
      celldata=["gcell"=>own_to_global(gids),"aggregate"=>own_aggregates])
    writevtk(Γ,"boudary")
    writevtk(Ω,"trian")
  end

  # No remotes

  np = (2,2)
  nc = (8,8)
  R = 0.4
  geo = disk(R,x0=Point(0.5,0.5))
  domain = (0, 1, 0, 1)
  bgmodel = CartesianDiscreteModel(ranks,np,domain,nc)
  cutgeo = cut(bgmodel,geo)
  strategy = AggregateCutCellsByThreshold(1.0)
  aggregates,aggregate_owner = distributed_aggregate(strategy,cutgeo)

  @test ! has_remote_aggregation(bgmodel,aggregates)

end

end # module
