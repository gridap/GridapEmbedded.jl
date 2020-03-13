
struct SubTriangulation{Dp,T}
  cell_to_points::Table{Int,Int32}
  cell_to_inoutcut::Vector{Int8}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dp,T}}
end

struct FacetSubTriangulation{Dp,T}
  facet_to_points::Table{Int,Int32}
  facet_to_normal::Vector{Point{Dp,T}}
  facet_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dp,T}}
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

function initial_sub_triangulation(grid::Grid,geom::AnalyticalGeometry)
  geom_x = discretize(geom,grid)
  initial_sub_triangulation(grid,geom_x)
end

function initial_sub_triangulation(grid::Grid,geom::DiscreteGeometry)
  ugrid = UnstructuredGrid(grid)
  _initial_sub_triangulation(ugrid,geom)
end

# Implementation of the Grid interface

# Helpers

struct MiniCell
  num_points::Int
  edge_to_points::Vector{Vector{Int}}
end

function _initial_sub_triangulation(grid::UnstructuredGrid,geom::DiscreteGeometry)

  cutgrid, cutgeom = _extract_grid_of_cut_cells(grid,geom)

  subtrian, subgeom = _simplexify_and_isolate_cells_in_cutgrid(cutgrid,cutgeom)

  subtrian, subgeom
end

function _extract_grid_of_cut_cells(grid,geom)

  p = _check_and_get_polytope(grid)
  table = LookupTable(p)
  cell_to_points = get_cell_nodes(grid)
  cell_to_inoutcut = _compute_in_out_or_cut(
    table,cell_to_points,geom.ls_to_point_to_value,geom.intersection)
  cutcell_to_cell = findall(cell_to_inoutcut .== CUT)
  cutgrid = GridPortion(grid,cutcell_to_cell)
  ls_to_point_to_value = [
    point_to_value[cutgrid.node_to_oldnode] for point_to_value in geom.ls_to_point_to_value ]
  cutgeom = DiscreteGeometry(geom.pmin,geom.pmax,geom.intersection,ls_to_point_to_value)

  cutgrid, cutgeom
end

function _check_and_get_polytope(grid)
  reffes = get_reffes(grid)
  @notimplementedif length(reffes) != 1
  reffe = first(reffes)
  order = 1
  @notimplementedif get_order(reffe) != order
  p = get_polytope(reffe)
  p
end

function _simplexify_and_isolate_cells_in_cutgrid(cutgrid,cutgeom)

  p = _check_and_get_polytope(cutgrid)

  ltcell_to_lpoints, simplex = simplexify(p)
  lpoint_to_lcoords = get_vertex_coordinates(p)
  _ensure_positive_jacobians!(ltcell_to_lpoints,lpoint_to_lcoords,simplex)

  out = _simplexify(
    get_node_coordinates(cutgrid),
    get_cell_nodes(cutgrid),
    ltcell_to_lpoints,
    lpoint_to_lcoords,
    cutgeom.ls_to_point_to_value,
    num_vertices(p),
    num_vertices(simplex))

  tcell_to_tpoints, tpoint_to_coords, tpoint_to_rcoords, ls_to_tpoint_to_value = out

  ntcells = length(tcell_to_tpoints)
  tcell_to_inoutcut = fill(Int8(CUT),ntcells)
  nltcells = length(ltcell_to_lpoints)
  tcell_to_cell = _setup_cell_to_bgcell(cutgrid.cell_to_oldcell,nltcells,ntcells)

  subtrian = SubTriangulation(
    tcell_to_tpoints,
    tcell_to_inoutcut,
    tcell_to_cell,
    tpoint_to_coords,
    tpoint_to_rcoords)

  subgeom = DiscreteGeometry(
    cutgeom.pmin,
    cutgeom.pmax,
    cutgeom.intersection,
    ls_to_tpoint_to_value)

  subtrian, subgeom
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

function _simplexify(
  point_to_coords,
  cell_to_points::Table,
  ltcell_to_lpoints,
  lpoint_to_lcoords,
  ls_to_point_to_value,
  nlpoints,
  nsp)

  ncells = length(cell_to_points)
  nltcells = length(ltcell_to_lpoints)
  ntcells = ncells*nltcells
  ntpoints = ncells*nlpoints

  tcell_to_tpoints_data = zeros(eltype(cell_to_points.data),nsp*ntcells)
  tcell_to_tpoints_ptrs = fill(eltype(cell_to_points.ptrs)(nsp),ntcells+1)
  length_to_ptrs!(tcell_to_tpoints_ptrs)
  tcell_to_tpoints = Table(tcell_to_tpoints_data,tcell_to_tpoints_ptrs)
  tpoint_to_coords = zeros(eltype(point_to_coords),ntpoints)
  tpoint_to_rcoords = zeros(eltype(point_to_coords),ntpoints)
  T = eltype(first(ls_to_point_to_value))
  ls_to_tpoint_to_value = [ zeros(T,ntpoints) for i in 1:length(ls_to_point_to_value)]

  tpoint = 0
  tcell = 0
  for cell in 1:ncells

    for ltcell in 1:nltcells
      tcell += 1
      q = tcell_to_tpoints.ptrs[tcell] - 1
      lpoints = ltcell_to_lpoints[ltcell]
      for (j,lpoint) in enumerate(lpoints)
        tcell_to_tpoints.data[q+j] = tpoint + lpoint
      end
    end

    a = cell_to_points.ptrs[cell]-1
    for lpoint in 1:nlpoints
      tpoint += 1
      point = cell_to_points.data[a+lpoint]
      coords = point_to_coords[point]
      rcoords = lpoint_to_lcoords[lpoint]
      tpoint_to_coords[tpoint] = coords
      tpoint_to_rcoords[tpoint] = rcoords
      for (i,point_to_val) in enumerate(ls_to_point_to_value)
        val = point_to_val[point]
        ls_to_tpoint_to_value[i][tpoint] = val
      end
    end

  end

  tcell_to_tpoints, tpoint_to_coords, tpoint_to_rcoords, ls_to_tpoint_to_value
end

@inline function _compute_in_out_or_cut(
  table::LookupTable,
  cell_to_points::Table,
  point_to_value::AbstractVector,
  cell::Integer)

  case = compute_case(cell_to_points,point_to_value,cell)
  table.case_to_inoutcut[case]
end

@inline function _compute_in_out_or_cut(
  table::LookupTable,
  cell_to_points::Table,
  ls_to_point_to_value::Vector{<:AbstractVector},
  cell::Integer,
  intersection::Bool=true)

  at_least_one_cut = false
  all_in = true
  all_out = true
  all_cut_or_in = true
  all_cut_or_out = true

  for point_to_value in ls_to_point_to_value
    inoutcut = _compute_in_out_or_cut(table,cell_to_points,point_to_value,cell)
    at_least_one_cut = at_least_one_cut || (inoutcut == CUT)
    all_in = all_in && (inoutcut == IN)
    all_out = all_out && (inoutcut == IN)
    all_cut_or_in = all_cut_or_in  && ( (inoutcut == CUT) || (inoutcut == IN) )
    all_cut_or_out = all_cut_or_out  && ( (inoutcut == CUT) || (inoutcut == OUT) )
  end

  if intersection
    if all_in
      r = IN
    elseif all_cut_or_in && at_least_one_cut
      r = CUT
    else
      r = OUT
    end
  else
    if all_out
      r = OUT
    elseif all_cut_or_out && at_least_one_cut
      r = CUT
    else
      r = IN
    end
  end

  r
end

function _compute_in_out_or_cut(
  table::LookupTable,
  cell_to_points::Table,
  point_to_value::AbstractVector)

  ncells = length(cell_to_points)
  cell_to_inoutcut = zeros(Int8,ncells)
  for cell in 1:ncells
    cell_to_inoutcut[cell] = _compute_in_out_or_cut(table,cell_to_points,point_to_value,cell)
  end
  cell_to_inoutcut
end

function _compute_in_out_or_cut(
  table::LookupTable,
  cell_to_points::Table,
  point_to_value::Vector{<:AbstractVector},
  intersection::Bool=true)

  ncells = length(cell_to_points)
  cell_to_inoutcut = zeros(Int8,ncells)
  for cell in 1:ncells
    cell_to_inoutcut[cell] = _compute_in_out_or_cut(table,cell_to_points,point_to_value,cell,intersection)
  end
  cell_to_inoutcut
end
