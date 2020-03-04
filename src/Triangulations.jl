
@inline function compute_in_out_or_cut(
  table::LookupTable,
  cell_to_points::Table,
  point_to_value::AbstractVector,
  cell::Integer)

  case = compute_case(cell_to_points,point_to_value,cell)
  table.case_to_inoutcut[case]
end

@inline function compute_in_out_or_cut(
  table::LookupTable,
  cell_to_points::Table,
  ls_to_point_to_value::Vector{<:AbstractVector},
  cell::Integer)

  for point_to_value in ls_to_point_to_value
    if OUT == compute_in_out_or_cut(table,cell_to_points,point_to_value,cell)
      return OUT
    end
  end
  for point_to_value in ls_to_point_to_value
    if CUT == compute_in_out_or_cut(table,cell_to_points,point_to_value,cell)
      return CUT
    end
  end
  return IN
end

function compute_in_out_or_cut(
  table::LookupTable,
  cell_to_points::Table,
  point_to_value::AbstractVector)

  ncells = length(cell_to_points)
  cell_to_inoutcut = zeros(Int8,ncells)
  for cell in 1:ncells
    cell_to_inoutcut[cell] = compute_in_out_or_cut(table,cell_to_points,point_to_value,cell)
  end
  cell_to_inoutcut
end

struct SubTriangulation{D,T}
  table::LookupTable{D,T}
  cell_to_points::Table{Int,Int32}
  cell_to_inoutcut::Vector{Int8}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{D,T}}
  ls_to_point_to_value::Vector{Vector{T}}
end

function initial_sub_triangulation(_grid::Grid,_ls_to_point_to_value::Vector{<:AbstractVector})
  grid = UnstructuredGrid(_grid)
  reffes = get_reffes(grid)
  @notimplementedif length(reffes) != 1
  reffe = first(reffes)
  order = 1
  @notimplementedif get_order(reffe) != order
  p = get_polytope(reffe)
  _table = LookupTable(p)
  _cell_to_points = get_cell_nodes(grid)
  _cell_to_inoutcut = compute_in_out_or_cut(_table,_cell_to_points,_ls_to_point_to_value)
  cutcell_to_cell = findall(_cell_to_inoutcut .== CUT)
  _cutgrid = GridPortion(grid,cutcell_to_cell)
  ls_to_point_to_value = [ point_to_value[_cutgrid.node_to_oldnode] for point_to_value in _ls_to_point_to_value ]
  cutgrid = simplexify(_cutgrid)
  ltcell_to_lpoints, simplex = simplexify(p)
  point_to_coords = get_node_coordinates(cutgrid)
  cell_to_points = get_cell_nodes(cutgrid)
  ncells = length(cell_to_points)
  cell_to_inoutcut = fill(Int8(CUT),ncells)
  nltcells = length(ltcell_to_lpoints)
  cell_to_bgcell = _setup_cell_to_bgcell(_cutgrid.cell_to_oldcell,nltcells,ncells)
  table = LookupTable(simplex)

  SubTriangulation(
    table,
    cell_to_points,
    cell_to_inoutcut,
    cell_to_bgcell,
    point_to_coords,
    ls_to_point_to_value)
end

function _setup_cell_to_bgcell(pcell_to_bgcell,nlcells,ncells)
  cell_to_bgcell = zeros(Int32,ncells)
  cell = 1
  for bgcell in pcell_to_bgcell
    for lcell in 1:nlcells
      cell_to_bgcell[cell] = bgcell
      cell += 1
    end
  end
  cell_to_bgcell
end

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
  nls = length(st.ls_to_point_to_value)
  nodaldata = [
    "ls_$(lpad(i,ceil(Int,log10(nls)),'0'))" => ls  for (i,ls) in enumerate(st.ls_to_point_to_value)  ]
  celldata = ["inoutcut"=>st.cell_to_inoutcut, "bgcell"=>st.cell_to_bgcell]
  write_vtk_file(ug,filename,
    celldata=celldata, nodaldata=nodaldata)
end
