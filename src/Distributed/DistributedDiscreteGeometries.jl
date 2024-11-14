struct DistributedDiscreteGeometry{A} <: GridapType
  geometries::A
end

local_views(a::DistributedDiscreteGeometry) = a.geometries

function DiscreteGeometry(φh::CellField,model::DistributedDiscreteModel;name::String="")
  geometries = map(local_views(φh),local_views(model)) do φh, model
    DiscreteGeometry(φh,model;name)
  end
  DistributedDiscreteGeometry(geometries)
end

function distributed_geometry(a::AbstractArray{<:DiscreteGeometry})
  DistributedDiscreteGeometry(a)
end

function discretize(a::AnalyticalGeometry,model::DistributedDiscreteModel)
  geometries = map(local_views(model)) do model
    discretize(a,model)
  end
  DistributedDiscreteGeometry(geometries)
end

function cut(
  cutter::Cutter,bgmodel::DistributedDiscreteModel,geom::DistributedDiscreteGeometry
)
  gids = get_cell_gids(bgmodel)
  cuts = map(local_views(bgmodel),local_views(geom)) do bgmodel,geom
    cut(cutter,bgmodel,geom)
  end
  @notimplementedif !isconsistent_bgcell_to_inoutcut(cuts,partition(gids))
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

function cut_facets(
  cutter::Cutter,bgmodel::DistributedDiscreteModel{Dc},geom::DistributedDiscreteGeometry
) where Dc
  gids = get_face_gids(bgmodel,Dc-1)
  cuts = map(local_views(bgmodel),local_views(geom)) do bgmodel,geom
    cut_facets(cutter,bgmodel,geom)
  end
  @notimplementedif !isconsistent_bgcell_to_inoutcut(cuts,partition(gids))
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

for TT in (:Triangulation,:SkeletonTriangulation,:BoundaryTriangulation,:EmbeddedBoundary)
  @eval begin
    function $TT(portion,cutgeo::DistributedEmbeddedDiscretization,cutinorout,geom::DistributedDiscreteGeometry)
      model = get_background_model(cutgeo)
      gids  = get_cell_gids(model)
      trians = map(local_views(cutgeo),local_views(geom),partition(gids)) do cutgeo, geom, gids
        $TT(portion,gids,cutgeo,cutinorout,geom)
      end
      DistributedTriangulation(trians,model)
    end
  end
end

function distributed_aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::DistributedEmbeddedDiscretization,
  geo::DistributedDiscreteGeometry,
  in_or_out = IN
)
  bgmodel = get_background_model(cut)
  facet_to_inoutcut = compute_bgfacet_to_inoutcut(bgmodel,geo)
  _distributed_aggregate_by_threshold(strategy.threshold,cut,geo,in_or_out,facet_to_inoutcut)
end

function compute_bgcell_to_inoutcut(
  cutgeo::DistributedEmbeddedDiscretization,geo::DistributedDiscreteGeometry
)
  map(local_views(cutgeo),local_views(geo)) do cutgeo, geo
    compute_bgcell_to_inoutcut(cutgeo,geo)
  end
end

function compute_bgfacet_to_inoutcut(
  cutter::Cutter, bgmodel::DistributedDiscreteModel, geo::DistributedDiscreteGeometry
)
  map(local_views(bgmodel),local_views(geo)) do model, geo
    compute_bgfacet_to_inoutcut(cutter,model,geo)
  end
end