
# TODO move to Gridap
function get_face_to_face(trian::BoundaryTriangulation)
  @abstractmethod
end

# TODO move to Gridap
function get_face_to_face(trian::GenericBoundaryTriangulation)
  trian.face_trian.cell_to_oldcell
end

struct EmbeddedFacetDiscretization{Dc,Dp,T} <: GridapType
  bgmodel::DiscreteModel{Dp,Dp}
  ls_to_subfacet_to_inout::Vector{Vector{Int8}}
  subfacets::SubTriangulation{Dc,Dp,T}
  ls_to_facet_to_inoutcut::Vector{Vector{Int8}}
  oid_to_ls::Dict{UInt,Int}
  geo::CSG.Geometry
end

function Gridap.BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization,geo::CSG.Geometry,in_or_out::Integer)

  bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(cut,geo)
  facet_to_mask = collect(Bool,bgfacet_to_inoutcut .== in_or_out)
  BoundaryTriangulation(cut.bgmodel,facet_to_mask)
end

function Gridap.BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization,geo::CSG.Geometry,in_or_out::CutInOrOut)

  bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(cut,geo)
  subfacet_to_inoutcut = reindex(bgfacet_to_inoutcut,cut.subfacets.cell_to_bgcell)


end

function compute_bgfacet_to_inoutcut(cut::EmbeddedFacetDiscretization,geo::CSG.Geometry)

  tree = get_tree(geo)

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    ls = cut.oid_to_ls[oid]
    cell_to_inoutcut = cut.ls_to_facet_to_inoutcut[ls]
    cell_to_inoutcut, name, meta
  end

  newtree = replace_data(identity,conversion,tree)
  compute_inoutcut(newtree)
end







struct BoundarySubTriangulationWrapper{Dc,Dp,T} <: Triangulation{Dc,Dp}
  facets::BoundaryTriangulation{Dc,Dp}
  subfacets::SubTriangulation{Dc,Dp,T}
  reffes::Vector{LagrangianRefFE{Dp}}
  cell_types::Vector{Int8}
  cell_ids
  cell_normals
  subfacet_to_facet_map

  function BoundarySubTriangulationWrapper(
    facets::BoundaryTriangulation{Dc,Dp},
    subfacets::SubTriangulation{Dc,Dp,T}) where {Dc,Dp,T}

    reffe = LagrangianRefFE(Float64,Simplex(Val{Dp}()),1)
    cell_types = fill(Int8(1),length(st.cell_to_points))
    reffes = [reffe]
    cell_ids = reindex(get_cell_id(facets),subfacets.cell_to_bgcell)
    cell_normals = reindex(get_normal_vector(facets),subfacets.cell_to_bgcell)
    subfacet_to_facet_map = _setup_subcell_to_cell_map(subfacets,reffe,cell_types)
    new{Dc,Dp,T}(facets,subfacets,reffes,cell_types,cell_ids,cell_normals)
  end
end

function get_node_coordinates(trian::BoundarySubTriangulationWrapper)
  trian.subfacets.point_to_coords
end

function get_cell_nodes(trian::BoundarySubTriangulationWrapper)
  trian.subfacets.cell_to_points
end

function get_cell_coordinates(trian::BoundarySubTriangulationWrapper)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_nodes(trian)
  LocalToGlobalArray(cell_to_nodes,node_to_coords)
end

function get_reffes(trian::BoundarySubTriangulationWrapper)
  trian.reffes
end

function get_cell_type(trian::BoundarySubTriangulationWrapper)
  trian.cell_types
end

function get_normal_vector(trian::BoundarySubTriangulationWrapper)
  trian.cell_normals
end

function get_cell_id(trian::BoundarySubTriangulationWrapper)
  trian.cell_ids
end

function restrict(f::AbstractArray,trian::BoundarySubTriangulationWrapper)
  compose_field_arrays(restrict(f,trian.facets), trian.subfacet_to_facet_map)
end


