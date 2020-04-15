
struct SubTriangulation{Dp,T} <: GridapType
  cell_to_points::Table{Int,Int32}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dp,T}}
end

function SubTriangulation(st::SubTriangulation,newcells::AbstractVector{<:Integer})
  cell_to_points = Table(reindex(st.cell_to_points,newcells))
  cell_to_bgcell = st.cell_to_bgcell[newcells]
  SubTriangulation(
    cell_to_points,
    cell_to_bgcell,
    st.point_to_coords,
    st.point_to_rcoords)
end


# Implementation of Triangulation interface

struct SubTriangulationWrapper{Dp,T} <: Triangulation{Dp,Dp}
  subcells::SubTriangulation{Dp,T}
  cell_types::Vector{Int8}
  reffes::Vector{LagrangianRefFE{Dp}}
  subcell_to_cell_map

  function SubTriangulationWrapper(st::SubTriangulation{Dp,T}) where {Dp,T}
    reffe = LagrangianRefFE(Float64,Simplex(Val{Dp}()),1)
    cell_types = fill(Int8(1),length(st.cell_to_points))
    reffes = [reffe]
    subcell_to_cell_map = _setup_subcell_to_cell_map(st,reffe,cell_types)
    new{Dp,T}(st,cell_types,reffes,subcell_to_cell_map)
  end
end

function _setup_subcell_to_cell_map(st,reffe,cell_types)
  subcell_to_rcoords = LocalToGlobalArray(st.cell_to_points,st.point_to_rcoords)
  subcell_to_shapefuns = CompressedArray([get_shapefuns(reffe)],cell_types)
  subcell_to_cell_map = lincomb(subcell_to_shapefuns,subcell_to_rcoords)
  subcell_to_cell_map
end

function get_node_coordinates(trian::SubTriangulationWrapper)
  trian.subcells.point_to_coords
end

function get_cell_nodes(trian::SubTriangulationWrapper)
  trian.subcells.cell_to_points
end

function get_reffes(trian::SubTriangulationWrapper)
  trian.reffes
end

function get_cell_type(trian::SubTriangulationWrapper)
  trian.cell_types
end

function get_cell_id(trian::SubTriangulationWrapper)
  trian.subcells.cell_to_bgcell
end

function get_cell_coordinates(trian::SubTriangulationWrapper)
  node_to_coords = get_node_coordinates(trian)
  cell_to_nodes = get_cell_nodes(trian)
  LocalToGlobalArray(cell_to_nodes,node_to_coords)
end

function restrict(f::AbstractArray, trian::SubTriangulationWrapper)
  compose_field_arrays(reindex(f,trian), trian.subcell_to_cell_map)
end

# API

function UnstructuredGrid(st::SubTriangulation{D}) where D
  reffe = LagrangianRefFE(Float64,Simplex(Val{D}()),1)
  cell_types = fill(Int8(1),length(st.cell_to_points))
  UnstructuredGrid(
    st.point_to_coords,
    st.cell_to_points,
    [reffe,],
    cell_types)
end

function writevtk(st::SubTriangulation,filename::String)
  ug = UnstructuredGrid(st)
  degree = 0
  quad = CellQuadrature(ug,degree)
  dV = integrate(1,ug,quad)
  write_vtk_file(ug,filename,celldata=[
    "bgcell"=>st.cell_to_bgcell,
    "dV"=>dV])
end

