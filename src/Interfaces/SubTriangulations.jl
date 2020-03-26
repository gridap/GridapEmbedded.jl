
struct SubTriangulation{Dp,T} <: GridapType
  cell_to_points::Table{Int,Int32}
  cell_to_inoutcut::Vector{Int8}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dp,T}}
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
    "inoutcut"=>st.cell_to_inoutcut,
    "bgcell"=>st.cell_to_bgcell,
    "dV"=>dV])
end

function split_in_out(st::SubTriangulation)
 st_in = take_in_or_out(st,IN)
 st_out = take_in_or_out(st,OUT)
 st_in, st_out 
end

function take_in_or_out(st::SubTriangulation{D},in_or_out) where D
  ntcells = 0
  for inoutcut in st.cell_to_inoutcut
    if  inoutcut == in_or_out
      ntcells += 1
    end
  end
  tcell_to_inoutcut = fill(Int8(in_or_out),ntcells)
  tcell_to_bgcell = zeros(Int32,ntcells)
  ntlpoints = D + 1
  tcell_to_points_ptrs = fill(Int32(ntlpoints),ntcells+1)
  length_to_ptrs!(tcell_to_points_ptrs)
  ndata = tcell_to_points_ptrs[end]-1
  tcell_to_points_data = zeros(Int,ndata)
  tcell_to_points = Table(tcell_to_points_data,tcell_to_points_ptrs)
  tcell = 0
  for (cell, inoutcut) in enumerate(st.cell_to_inoutcut)
    if  inoutcut == in_or_out
      tcell += 1
      tcell_to_bgcell[tcell] = st.cell_to_bgcell[cell]
      a = tcell_to_points.ptrs[tcell]-1
      b = st.cell_to_points.ptrs[cell]-1
      for i in 1:ntlpoints
        tcell_to_points.data[a+i] = st.cell_to_points.data[b+i]
      end
    end
  end
  SubTriangulation(
    tcell_to_points,
    tcell_to_inoutcut,
    tcell_to_bgcell,
    st.point_to_coords,
    st.point_to_rcoords)
end

