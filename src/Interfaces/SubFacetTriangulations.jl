
struct SubFacetData{Dp,T,Tn} <: GridapType
  facet_to_points::Table{Int32,Vector{Int32},Vector{Int32}}
  facet_to_normal::Vector{Point{Dp,Tn}}
  facet_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dp,T}}
end

function SubFacetData(
  st::SubFacetData,newfacets::AbstractVector{<:Integer},orientation::AbstractVector)
  facet_to_points = Table(lazy_map(Reindex(st.facet_to_points),newfacets))
  facet_to_normal = st.facet_to_normal[newfacets] .* orientation
  facet_to_bgcell = st.facet_to_bgcell[newfacets]
  SubFacetData(
    facet_to_points,
    facet_to_normal,
    facet_to_bgcell,
    st.point_to_coords,
    st.point_to_rcoords)
end

# Implementation of the Gridap.Triangulation interface

struct SubFacetTriangulation{Dc,Dp,T,Tn} <: Grid{Dc,Dp}
  subfacets::SubFacetData{Dp,T,Tn}
  bgtrian::Triangulation
  facet_types::Vector{Int8}
  reffes::Vector{LagrangianRefFE{Dc}}
  facet_ref_map

  function SubFacetTriangulation(st::SubFacetData{Dp,T,Tn},bgtrian::Triangulation) where {Dp,T,Tn}
    Dc = Dp - 1
    reffe = LagrangianRefFE(Float64,Simplex(Val{Dc}()),1)
    facet_types = fill(Int8(1),length(st.facet_to_points))
    reffes = [reffe]
    face_ref_map = _setup_facet_ref_map(st,reffe,facet_types)
    new{Dc,Dp,T,Tn}(st,bgtrian,facet_types,reffes,face_ref_map)
  end
end

function _setup_facet_ref_map(st,reffe,facet_types)
  facet_to_points = st.facet_to_points
  point_to_rcoords = st.point_to_rcoords
  facet_to_rcoords = lazy_map(Broadcasting(Reindex(point_to_rcoords)),facet_to_points)
  facet_to_shapefuns = expand_cell_data([get_shapefuns(reffe)],facet_types)
  facet_to_ref_map = lazy_map(linear_combination,facet_to_rcoords,facet_to_shapefuns)
  facet_to_ref_map
end

function compress_contributions(cell_mat,trian::SubFacetTriangulation)
  cell_to_bgcell = get_cell_to_bgcell(trian)
  ccell_mat = compress_contributions(cell_mat,cell_to_bgcell)
  ccell_mat
end

function compress_ids(cell_ids,trian::SubFacetTriangulation)
  cell_to_bgcell = get_cell_to_bgcell(trian)
  compress_ids(cell_ids,cell_to_bgcell)
end

# Triangulation API

get_node_coordinates(trian::SubFacetTriangulation) = trian.subfacets.point_to_coords
get_cell_node_ids(trian::SubFacetTriangulation) = trian.subfacets.facet_to_points
get_reffes(trian::SubFacetTriangulation) = trian.reffes
get_cell_type(trian::SubFacetTriangulation) = trian.facet_types
TriangulationStyle(::Type{<:SubFacetTriangulation}) = SubTriangulation()
get_background_triangulation(trian::SubFacetTriangulation) = trian.bgtrian
get_cell_to_bgcell(trian::SubFacetTriangulation) = trian.subfacets.facet_to_bgcell
get_cell_ref_map(trian::SubFacetTriangulation) = trian.facet_ref_map
get_facet_normal(trian::SubFacetTriangulation) = lazy_map(ConstantField,trian.subfacets.facet_to_normal)

# API

function UnstructuredGrid(st::SubFacetData{Dp}) where Dp
  Dc = Dp -1
  reffe = LagrangianRefFE(Float64,Simplex(Val{Dc}()),1)
  cell_types = fill(Int8(1),length(st.facet_to_points))
  UnstructuredGrid(
    st.point_to_coords,
    st.facet_to_points,
    [reffe,],
    cell_types)
end

function Visualization.visualization_data(st::SubFacetData,filename::String;celldata=Dict())
  ug = UnstructuredGrid(st)
  degree = 0
  quad = CellQuadrature(ug,degree)
  dS = integrate(1,quad)
  newcelldata = ["bgcell"=>st.facet_to_bgcell,"dS"=>dS,"normal"=>st.facet_to_normal]
  _celldata = Dict()
  for (k,v) in celldata
    _celldata[k] = v
  end
  for (k,v) in newcelldata
    _celldata[k] = v
  end
  (Visualization.VisualizationData(ug,filename,celldata=_celldata),)
end

function merge_sub_face_data(ls_to_subfacets,ls_to_ls_to_facet_to_inout)

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

function Base.empty(st::SubFacetData{Dp,T,Tn}) where {Dp,T,Tn}

  facet_to_points = Table(Int32[],Int32[1,])
  facet_to_normal = Point{Dp,Tn}[]
  facet_to_bgcell = Int32[]
  point_to_coords = Point{Dp,T}[]
  point_to_rcoords = Point{Dp,T}[]

  SubFacetData(
    facet_to_points,
    facet_to_normal,
    facet_to_bgcell,
    point_to_coords,
    point_to_rcoords)
end

function Base.append!(a::SubFacetData{D},b::SubFacetData{D}) where D

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
