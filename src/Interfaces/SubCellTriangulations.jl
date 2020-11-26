
struct SubCellData{Dr,Dp,T} <: GridapType
  cell_to_points::Table{Int32,Vector{Int32},Vector{Int32}}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dr,T}}
end

function SubCellData(st::SubCellData,newcells::AbstractVector{<:Integer})
  cell_to_points = Table(lazy_map(Reindex(st.cell_to_points),newcells))
  cell_to_bgcell = st.cell_to_bgcell[newcells]
  SubCellData(
    cell_to_points,
    cell_to_bgcell,
    st.point_to_coords,
    st.point_to_rcoords)
end

# Implementation of Triangulation interface

struct SubCellTriangulation{Dp,T} <: Grid{Dp,Dp}
  subcells::SubCellData{Dp,Dp,T}
  bgtrian::Triangulation
  cell_types::Vector{Int8}
  reffes::Vector{LagrangianRefFE{Dp}}
  cell_ref_map::AbstractArray{<:Field}

  function SubCellTriangulation(st::SubCellData{Dp,Dp,T},bgtrian::Triangulation) where {Dp,T}
    reffe = LagrangianRefFE(Float64,Simplex(Val{Dp}()),1)
    cell_types = fill(Int8(1),length(st.cell_to_points))
    reffes = [reffe]
    cell_to_ref_map = _setup_cell_ref_map(st,reffe,cell_types)
    new{Dp,T}(st,bgtrian,cell_types,reffes,cell_to_ref_map)
  end
end

function _setup_cell_ref_map(st,reffe,cell_types)
  cell_to_points = st.cell_to_points
  point_to_rcoords = st.point_to_rcoords
  cell_to_rcoords = lazy_map(Broadcasting(Reindex(point_to_rcoords)),cell_to_points)
  cell_to_shapefuns = expand_cell_data([get_shapefuns(reffe)],cell_types)
  cell_to_ref_map = lazy_map(linear_combination,cell_to_rcoords,cell_to_shapefuns)
  cell_to_ref_map
end

# Triangulation API

get_node_coordinates(trian::SubCellTriangulation) = trian.subcells.point_to_coords
get_cell_nodes(trian::SubCellTriangulation) = trian.subcells.cell_to_points
get_reffes(trian::SubCellTriangulation) = trian.reffes
get_cell_type(trian::SubCellTriangulation) = trian.cell_types
TriangulationStyle(::Type{<:SubCellTriangulation}) = SubTriangulation()
get_background_triangulation(trian::SubCellTriangulation) = trian.bgtrian
get_cell_id(trian::SubCellTriangulation) = trian.subcells.cell_to_bgcell
get_cell_ref_map(trian::SubCellTriangulation) = trian.cell_ref_map

# API

function UnstructuredGrid(st::SubCellData{D}) where D
  reffe = LagrangianRefFE(Float64,Simplex(Val{D}()),1)
  cell_types = fill(Int8(1),length(st.cell_to_points))
  UnstructuredGrid(
    st.point_to_coords,
    st.cell_to_points,
    [reffe,],
    cell_types)
end

function Visualization.visualization_data(st::SubCellData,filename::String;celldata=Dict())
  ug = UnstructuredGrid(st)
  degree = 0
  quad = CellQuadrature(ug,degree)
  dV = integrate(1,quad)
  newcelldata = ["bgcell"=>st.cell_to_bgcell, "dV"=>dV]
  _celldata = Dict()
  for (k,v) in celldata
    _celldata[k] = v
  end
  for (k,v) in newcelldata
    _celldata[k] = v
  end
  (Visualization.VisualizationData(ug,filename,celldata=_celldata),)
end

