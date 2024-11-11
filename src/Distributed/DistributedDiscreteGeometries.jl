struct DistributedDiscreteGeometry{A} <: GridapType
  geometries::A
end

local_views(a::DistributedDiscreteGeometry) = a.geometries

function _get_values_at_owned_coords(φh,model::DistributedDiscreteModel{Dc,Dp}) where {Dc,Dp}
  @assert DomainStyle(φh) == ReferenceDomain()
  gids = get_cell_gids(model)
  values = map(local_views(φh),local_views(model),local_views(gids)) do φh, model, gids
    # Maps from the no-ghost model to the original model
    own_model = remove_ghost_cells(model,gids)
    own_to_local_node = get_face_to_parent_face(own_model,0)
    local_to_own_node = find_inverse_index_map(own_to_local_node,num_nodes(model))
    own_to_local_cell = get_face_to_parent_face(own_model,Dc)

    # Cell-to-node map for the original model
    # topo = get_grid_topology(model)
    # c2n_map = get_faces(topo,Dc,0)
    c2n_map = collect1d(get_cell_node_ids(model))

    # Cell-wise node coordinates (in ReferenceDomain coordinates)
    cell_reffe = get_cell_reffe(model)
    cell_node_coords = lazy_map(get_node_coordinates,cell_reffe)

    φh_data = CellData.get_data(φh)
    T = return_type(testitem(CellData.get_data(φh)),testitem(testitem(cell_node_coords)))
    values  = Vector{T}(undef,num_nodes(own_model))
    touched = fill(false,num_nodes(model))

    cell_node_coords_cache = array_cache(cell_node_coords)
    for cell in own_to_local_cell # For each owned cell
      field = φh_data[cell]
      node_coords = getindex!(cell_node_coords_cache,cell_node_coords,cell)
      for (iN,node) in enumerate(c2n_map[cell]) # Go over local nodes
        own_node = local_to_own_node[node]
        if (own_node != 0) && !touched[node] # Compute value if suitable
          values[own_node] = field(node_coords[iN])
          touched[node] = true
        end
      end
    end
    return values
  end
  return values
end

function DiscreteGeometry(φh::CellField,model::DistributedDiscreteModel;name::String="")
  φ_values = _get_values_at_owned_coords(φh,model)
  gids = get_cell_gids(model)
  geometries = map(local_views(model),local_views(gids),local_views(φ_values)) do model,gids,loc_φ
    ownmodel = remove_ghost_cells(model,gids)
    point_to_coords = collect1d(get_node_coordinates(ownmodel))
    DiscreteGeometry(loc_φ,point_to_coords;name)
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