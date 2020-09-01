
struct FacetSubTriangulation{Dp,T} <: GridapType
  facet_to_points::Table{Int,Vector{Int},Vector{Int32}}
  facet_to_normal::Vector{Point{Dp,T}}
  facet_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dp,T}}
end

function FacetSubTriangulation(
  st::FacetSubTriangulation,newfacets::AbstractVector{<:Integer},orientation::AbstractVector)
  facet_to_points = Table(reindex(st.facet_to_points,newfacets))
  facet_to_normal = st.facet_to_normal[newfacets] .* orientation
  facet_to_bgcell = st.facet_to_bgcell[newfacets]
  FacetSubTriangulation(
    facet_to_points,
    facet_to_normal,
    facet_to_bgcell,
    st.point_to_coords,
    st.point_to_rcoords)
end

# Implementation of the Gridap.Triangulation interface

struct FacetSubTriangulationWrapper{Dc,Dp,T} <: Triangulation{Dc,Dp}
  subfacets::FacetSubTriangulation{Dp,T}
  cell_types::Vector{Int8}
  reffes::Vector{LagrangianRefFE{Dc}}
  face_to_cell_map

  function FacetSubTriangulationWrapper(st::FacetSubTriangulation{Dp,T}) where {Dp,T}
    Dc = Dp - 1
    reffe = LagrangianRefFE(Float64,Simplex(Val{Dc}()),1)
    cell_types = fill(Int8(1),length(st.facet_to_points))
    reffes = [reffe]
    face_to_cell_map = _setup_face_to_cell_map(st,reffe,cell_types)
    new{Dc,Dp,T}(st,cell_types,reffes,face_to_cell_map)
  end
end

function _setup_face_to_cell_map(st,reffe,cell_types)
  facet_to_rcoords = LocalToGlobalArray(st.facet_to_points,st.point_to_rcoords)
  facet_to_shapefuns = CompressedArray([get_shapefuns(reffe)],cell_types)
  face_to_cell_map = lincomb(facet_to_shapefuns,facet_to_rcoords)
  face_to_cell_map
end

function get_node_coordinates(trian::FacetSubTriangulationWrapper)
  trian.subfacets.point_to_coords
end

function get_cell_nodes(trian::FacetSubTriangulationWrapper)
  trian.subfacets.facet_to_points
end

function get_reffes(trian::FacetSubTriangulationWrapper)
  trian.reffes
end

function get_cell_type(trian::FacetSubTriangulationWrapper)
  trian.cell_types
end

function get_normal_vector(trian::FacetSubTriangulationWrapper)
  cell_map = get_cell_map(trian)
  a = trian.subfacets.facet_to_normal
  GenericCellField(a,cell_map)
end

function get_face_to_cell(trian::FacetSubTriangulationWrapper)
  trian.subfacets.facet_to_bgcell
end

function get_face_to_cell_map(trian::FacetSubTriangulationWrapper)
  trian.face_to_cell_map
end

function get_cell_coordinates(trian::FacetSubTriangulationWrapper)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_nodes(trian)
  LocalToGlobalArray(cell_to_nodes,node_to_coords)
end

function restrict(f::AbstractArray, trian::FacetSubTriangulationWrapper)
  compose_field_arrays(reindex(f,trian), get_face_to_cell_map(trian))
end

function get_cell_id(trian::FacetSubTriangulationWrapper)
  get_face_to_cell(trian)
end

# API

function UnstructuredGrid(st::FacetSubTriangulation{Dp}) where Dp
  Dc = Dp -1
  reffe = LagrangianRefFE(Float64,Simplex(Val{Dc}()),1)
  cell_types = fill(Int8(1),length(st.facet_to_points))
  UnstructuredGrid(
    st.point_to_coords,
    st.facet_to_points,
    [reffe,],
    cell_types)
end

function writevtk(st::FacetSubTriangulation,filename::String,celldata=[])
  ug = UnstructuredGrid(st)
  degree = 0
  quad = CellQuadrature(ug,degree)
  dS = integrate(1,ug,quad)

  newcelldata = [
    "normal"=>st.facet_to_normal,
    "bgcell"=>st.facet_to_bgcell,
    "dS"=>dS]

  _celldata = vcat(celldata,newcelldata)

  write_vtk_file(ug,filename,celldata=_celldata)
end

function merge_facet_sub_triangulations(ls_to_subfacets,ls_to_ls_to_facet_to_inout)

  fst = empty(first(ls_to_subfacets))
  ls_to_facet_to_inout = [ empty(first(ls_to_f_to_i))  for ls_to_f_to_i in ls_to_ls_to_facet_to_inout  ]
  for (i,fst_i) in enumerate(ls_to_subfacets)
    append!(fst,fst_i)
    for (j,facet_to_inout) in enumerate(ls_to_ls_to_facet_to_inout[i])
      append!(ls_to_facet_to_inout[j],facet_to_inout)
    end
  end
  fst, ls_to_facet_to_inout
end

function Base.empty(st::FacetSubTriangulation{Dp,T}) where {Dp,T}

  facet_to_points = Table(Int[],Int32[1,])
  facet_to_normal = Point{Dp,T}[]
  facet_to_bgcell = Int32[]
  point_to_coords = Point{Dp,T}[]
  point_to_rcoords = Point{Dp,T}[]

  FacetSubTriangulation(
    facet_to_points,
    facet_to_normal,
    facet_to_bgcell,
    point_to_coords,
    point_to_rcoords)
end

function Base.append!(a::FacetSubTriangulation{D},b::FacetSubTriangulation{D}) where D

  o = length(a.point_to_coords)

  append!(a.facet_to_normal, b.facet_to_normal)
  append!(a.facet_to_bgcell, b.facet_to_bgcell)
  append!(a.point_to_coords, b.point_to_coords)
  append!(a.point_to_rcoords, b.point_to_rcoords)

  nini = length(a.facet_to_points.data)+1
  append!(a.facet_to_points.data, b.facet_to_points.data)
  nend = length(a.facet_to_points.data)
  append_ptrs!(a.facet_to_points.ptrs,b.facet_to_points.ptrs)
  for i in nini:nend
    a.facet_to_points.data[i] += o
  end

  a
end

