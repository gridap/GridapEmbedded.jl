
macro publish(mod,name)
  quote
    using GridapEmbedded.$mod: $name; export $name
  end
end

@publish CSG get_metadata
@publish CSG get_geometry
@publish CSG get_geometry_names

@publish Interfaces cut
@publish Interfaces cut_facets
@publish Interfaces compute_bgcell_to_inoutcut
@publish Interfaces compute_bgfacet_to_inoutcut
@publish Interfaces EmbeddedDiscretization
@publish Interfaces EmbeddedBoundary
@publish Interfaces GhostSkeleton
@publish Interfaces IN
@publish Interfaces OUT

@publish LevelSetCutters LevelSetCutter
@publish LevelSetCutters AnalyticalGeometry
@publish LevelSetCutters popcorn
@publish LevelSetCutters doughnut
@publish LevelSetCutters tube
@publish LevelSetCutters olympic_rings
@publish LevelSetCutters sphere
@publish LevelSetCutters disk
@publish LevelSetCutters cylinder
@publish LevelSetCutters plane
@publish LevelSetCutters cube
@publish LevelSetCutters square
@publish LevelSetCutters quadrilateral

@publish AgFEM AgFEMSpace
@publish AgFEM aggregate
@publish AgFEM color_aggregates
@publish AgFEM AggregateCutCellsByThreshold
@publish AgFEM AggregateAllCutCells
