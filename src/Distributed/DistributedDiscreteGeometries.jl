struct DistributedDiscreteGeometry{A} <: GridapType
  geometries::A
end

local_views(a::DistributedDiscreteGeometry) = a.geometries

function DistributedDiscreteGeometry(φh::CellField,model::DistributedDiscreteModel)
  gids = get_cell_gids(model)
  geometries = map(local_views(model),local_views(gids),local_views(φh)) do model,gids,φh
    ownmodel = remove_ghost_cells(model,gids)
    point_to_coords = collect1d(get_node_coordinates(ownmodel))
    DiscreteGeometry(φh(point_to_coords),point_to_coords)
  end
  DistributedDiscreteGeometry(geometries)
end

function get_geometry(a::AbstractArray{<:DiscreteGeometry})
  DistributedDiscreteGeometry(a)
end

function discretize(a::AnalyticalGeometry,model::DistributedDiscreteModel)
  gids = get_cell_gids(model)
  geometries = map(local_views(model),local_views(gids)) do model,gids
    ownmodel = remove_ghost_cells(model,gids)
    point_to_coords = collect1d(get_node_coordinates(ownmodel))
    discretize(a,point_to_coords)
  end
  DistributedDiscreteGeometry(geometries)
end

function cut(cutter::Cutter,bgmodel::DistributedDiscreteModel,geom::DistributedDiscreteGeometry)
  gids = get_cell_gids(bgmodel)
  cuts = map(local_views(bgmodel),local_views(gids),local_views(geom)) do bgmodel,gids,geom
    ownmodel = remove_ghost_cells(bgmodel,gids)
    cutgeo = cut(cutter,ownmodel,geom)
    change_bgmodel(cutgeo,bgmodel,own_to_local(gids))
  end
  consistent_bgcell_to_inoutcut!(cuts,gids)
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

function cut_facets(cutter::Cutter,bgmodel::DistributedDiscreteModel,geom::DistributedDiscreteGeometry)
  D = map(num_dims,local_views(bgmodel)) |> PartitionedArrays.getany
  cell_gids = get_cell_gids(bgmodel)
  facet_gids = get_face_gids(bgmodel,D-1)
  cuts = map(
    local_views(bgmodel),
    local_views(cell_gids),
    local_views(facet_gids),
    local_views(geom)) do bgmodel,cell_gids,facet_gids,geom
    ownmodel = remove_ghost_cells(bgmodel,cell_gids)
    facet_to_pfacet = get_face_to_parent_face(ownmodel,D-1)
    cutfacets = cut_facets(cutter,ownmodel,geom)
    cutfacets = change_bgmodel(cutfacets,bgmodel,facet_to_pfacet)
    remove_ghost_subfacets(cutfacets,facet_gids)
  end
  consistent_bgfacet_to_inoutcut!(cuts,facet_gids)
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

function distributed_embedded_triangulation(
  T,
  cutgeo::DistributedEmbeddedDiscretization,
  cutinorout,
  geom::DistributedDiscreteGeometry)

  trians = map(local_views(cutgeo),local_views(geom)) do lcutgeo,lgeom
    T(lcutgeo,cutinorout,lgeom)
  end
  bgmodel = get_background_model(cutgeo)
  DistributedTriangulation(trians,bgmodel)
end

function distributed_aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::DistributedEmbeddedDiscretization,
  geo::DistributedDiscreteGeometry,
  in_or_out=IN)

  bgmodel = get_background_model(cut)
  facet_to_inoutcut = compute_bgfacet_to_inoutcut(bgmodel,geo)
  _distributed_aggregate_by_threshold(strategy.threshold,cut,geo,in_or_out,facet_to_inoutcut)
end

function compute_bgcell_to_inoutcut(cutgeo::DistributedEmbeddedDiscretization,geo::DistributedDiscreteGeometry)
  map(local_views(cutgeo),local_views(geo)) do cutgeo,geo
    compute_bgcell_to_inoutcut(cutgeo,geo)
  end
end

function compute_bgfacet_to_inoutcut(
  cutter::Cutter,
  bgmodel::DistributedDiscreteModel,
  geo::DistributedDiscreteGeometry)

  gids = get_cell_gids(bgmodel)
  bgf_to_ioc = map(local_views(bgmodel),local_views(gids),local_views(geo)) do model,gids,geo
    ownmodel = remove_ghost_cells(model,gids)
    compute_bgfacet_to_inoutcut(cutter,ownmodel,geo)
  end
  compute_bgfacet_to_inoutcut(bgmodel,bgf_to_ioc)
end