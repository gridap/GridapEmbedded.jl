
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
  ls_to_facet_to_inoutcut::Vector{Vector{Int8}}
  subfacets::SubTriangulation{Dc,Dp,T}
  ls_to_subfacet_to_inout::Vector{Vector{Int8}}
  oid_to_ls::Dict{UInt,Int}
  geo::CSG.Geometry
end

function BoundaryTriangulation(cut::EmbeddedFacetDiscretization,tags,geo::CSG.Geometry,in_or_out)
  facets = BoundaryTriangulation(cut.bgmodel,tags)
  BoundaryTriangulation(cut,facets,geo,in_or_out)
end

function BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization,facets::BoundaryTriangulation,geo::CSG.Geometry,in_or_out::Tuple)

  trian1 = BoundaryTriangulation(cut,facets,geo,in_or_out[1])
  trian2 = BoundaryTriangulation(cut,facets,geo,in_or_out[2])
  lazy_append(trian1,trian2)
end

function BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization,facets::BoundaryTriangulation,geo::CSG.Geometry,in_or_out::Integer)

  facet_to_bgfacet = get_face_to_face(facets)
  bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(cut,geo)
  bgfacet_to_mask = collect(Bool,bgfacet_to_inoutcut .== in_or_out)
  facet_to_mask = reindex(bgfacet_to_mask,facet_to_bgfacet)

  TriangulationPortion(facets,findall(facet_to_mask))
end

function BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization,facets::BoundaryTriangulation,geo::CSG.Geometry,in_or_out::CutInOrOut)

  facet_to_bgfacet = get_face_to_face(facets)
  n_bgfacets = num_facets(cut.bgmodel)
  bgfacet_to_facet = zeros(Int,n_bgfacets)
  bgfacet_to_facet[facet_to_bgfacet] .= 1:length(facet_to_bgfacet)

  bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(cut,geo)
  subfacet_to_inoutcut = reindex(bgfacet_to_inoutcut,cut.subfacets.cell_to_bgcell)
  _subfacet_to_facet = reindex(bgfacet_to_facet,cut.subfacets.cell_to_bgcell)

  subfacet_to_inout = compute_subfacet_to_inout(cut,geo)
  pred(a,b,c) = c != 0 && a==CUT && b==in_or_out.in_or_out
  mask = apply( pred, subfacet_to_inoutcut, subfacet_to_inout, _subfacet_to_facet )
  newsubfacets = findall(mask)
  subfacets = SubTriangulation(cut.subfacets,newsubfacets)
  subfacet_to_facet = bgfacet_to_facet[subfacets.cell_to_bgcell]

  BoundarySubTriangulationWrapper(facets,subfacets,subfacet_to_facet)
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

function compute_subfacet_to_inout(cut::EmbeddedFacetDiscretization,geo::CSG.Geometry)

  tree = get_tree(geo)

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    ls = cut.oid_to_ls[oid]
    cell_to_inoutcut = cut.ls_to_subfacet_to_inout[ls]
    cell_to_inoutcut, name, meta
  end

  newtree = replace_data(identity,conversion,tree)
  compute_inoutcut(newtree)
end

struct BoundarySubTriangulationWrapper{Dc,Dp,T} <: Triangulation{Dc,Dp}
  facets::BoundaryTriangulation{Dc,Dp}
  subfacets::SubTriangulation{Dc,Dp,T}
  subfacet_to_facet::AbstractArray
  reffes::Vector{LagrangianRefFE{Dc}}
  cell_types::Vector{Int8}
  cell_ids
  cell_normals
  subfacet_to_facet_map

  function BoundarySubTriangulationWrapper(
    facets::BoundaryTriangulation{Dc,Dp},
    subfacets::SubTriangulation{Dc,Dp,T},
    subfacet_to_facet::AbstractArray) where {Dc,Dp,T}

    reffe = LagrangianRefFE(Float64,Simplex(Val{Dc}()),1)
    cell_types = fill(Int8(1),length(subfacets.cell_to_points))
    reffes = [reffe]
    cell_ids = reindex(get_cell_id(facets),subfacet_to_facet)
    cell_normals = reindex(get_normal_vector(facets),subfacet_to_facet)
    subfacet_to_facet_map = _setup_subcell_to_cell_map(subfacets,reffe,cell_types)
    new{Dc,Dp,T}(facets,subfacets,subfacet_to_facet,reffes,cell_types,cell_ids,cell_normals)
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
  g = restrict(f,trian.facets)
  compose_field_arrays(reindex(g,trian.subfacet_to_facet),trian.subfacet_to_facet_map)
end

