function cut_sub_triangulation(st::SubTriangulation{Dc,T},ls_to_point_to_value) where {Dc,T}
  _st = st

  _ls_to_point_to_value = [ point_to_value for point_to_value in ls_to_point_to_value ]
  ls_to_cell_to_inout = Vector{Int8}[]
  
  while length(_ls_to_point_to_value)>0
    @show "i"
    point_to_value = popfirst!(_ls_to_point_to_value)
    out = _cut_sub_triangulation(_st,point_to_value,_ls_to_point_to_value,ls_to_cell_to_inout)
    _st, _ls_to_point_to_value, ls_to_cell_to_inout, cell_to_inout = out
    push!(ls_to_cell_to_inout,cell_to_inout)
  end

  _st, ls_to_cell_to_inout
end

function _cut_sub_triangulation(st::SubTriangulation,point_to_value,i_to_point_to_value,j_to_cell_to_inout)

  refcell, reffacet, cell_table, facet_table = _setup_tables(st)
  nrcells, nrpoints = _cut_sub_triangulation_count(st,cell_table,point_to_value)
  rst, i_to_rpoint_to_value, j_to_rcell_to_inout, rcell_to_inout = _allocate_new_sub_triangulation(
    st,refcell,nrcells,nrpoints,i_to_point_to_value,j_to_cell_to_inout)
  _fill_new_sub_triangulation!(
    rst,
    i_to_rpoint_to_value,
    j_to_rcell_to_inout,
    rcell_to_inout,
    cell_table,
    refcell,
    st,
    point_to_value,
    i_to_point_to_value,
    j_to_cell_to_inout)
  rst, i_to_rpoint_to_value, j_to_rcell_to_inout, rcell_to_inout
end

struct MiniCell
  num_points::Int
  edge_to_points::Vector{Vector{Int}}
end

function _setup_tables(st::SubTriangulation{D}) where D
  p = Simplex(Val{D}())
  cell = MiniCell(num_vertices(p),get_faces(p,1,0))
  cell_table = LookupTable(p)
  pf = Simplex(Val{D-1}())
  facet = MiniCell(num_vertices(pf),get_faces(pf,1,0))
  facet_table = LookupTable(pf)
  cell, facet, cell_table, facet_table
end

function _cut_sub_triangulation_count(st::SubTriangulation,cell_table,point_to_value)
  nrcells = 0
  nrpoints = 0
  for cell in 1:length(st.cell_to_points)
    case = compute_case(st.cell_to_points,point_to_value,cell)
    nrcells += length(cell_table.case_to_subcell_to_inout[case])
    nrpoints += length(cell_table.case_to_point_to_coordinates[case])
  end
  nrcells, nrpoints
end

function _allocate_new_sub_triangulation(
  st::SubTriangulation{Dc,T},refcell,nrcells,nrpoints,i_to_point_to_value,j_to_cell_to_inout) where {Dc,T}

  nlp = refcell.num_points
  rcell_to_rpoints_data = zeros(eltype(st.cell_to_points.data),nlp*nrcells)
  rcell_to_rpoints_ptrs = fill(eltype(st.cell_to_points.ptrs)(nlp),nrcells+1)
  length_to_ptrs!(rcell_to_rpoints_ptrs)
  rcell_to_rpoints = Table(rcell_to_rpoints_data,rcell_to_rpoints_ptrs)
  rcell_to_bgcell = zeros(eltype(st.cell_to_bgcell),nrcells)
  rpoint_to_coords = zeros(Point{Dc,T},nrpoints)
  rpoint_to_rcoords = zeros(Point{Dc,T},nrpoints)
  i_to_rpoint_to_value = [zeros(T,nrpoints) for i in 1:length(i_to_point_to_value)]
  j_to_rcell_to_inout = [zeros(eltype(k),nrcells) for k in j_to_cell_to_inout]
  rcell_to_inout = zeros(eltype(eltype(j_to_rcell_to_inout)),nrcells)

  _st = SubTriangulation(
    rcell_to_rpoints,
    rcell_to_bgcell,
    rpoint_to_coords,
    rpoint_to_rcoords)

  _st, i_to_rpoint_to_value, j_to_rcell_to_inout, rcell_to_inout
end

function _fill_new_sub_triangulation!(
  rst::SubTriangulation,
  i_to_rpoint_to_value,
  j_to_rcell_to_inout,
  rcell_to_inout,
  cell_table,
  refcell,
  st,
  point_to_value,
  i_to_point_to_value,
  j_to_cell_to_inout)

  rcell = 0
  rpoint = 0
  rfacet = 0
  q = 0
  z = 0
  for cell in 1:length(st.cell_to_points)

    pointoffset = rpoint
    a = st.cell_to_points.ptrs[cell]-1
    for lpoint in 1:refcell.num_points
      rpoint += 1
      point = st.cell_to_points.data[a+lpoint]
      coords = st.point_to_coords[point]
      rst.point_to_coords[rpoint] = coords
      rcoords = st.point_to_rcoords[point]
      rst.point_to_rcoords[rpoint] = rcoords
      for (ils, point_to_val) in enumerate(i_to_point_to_value)
        val = point_to_val[point]
        i_to_rpoint_to_value[ils][rpoint] = val
      end
    end

    case = compute_case(st.cell_to_points,point_to_value,cell)

    if CUT == cell_table.case_to_inoutcut[case]
      for (ledge, lpoints) in enumerate(refcell.edge_to_points)
        point1 = st.cell_to_points.data[a+lpoints[1]]
        point2 = st.cell_to_points.data[a+lpoints[2]]
        v1 = point_to_value[point1]
        v2 = point_to_value[point2]
        if isout(v1) != isout(v2)
          rpoint += 1
          w1 = abs(v1)
          w2 = abs(v2)
          c1 = w1/(w1+w2)
          p1 = st.point_to_coords[point1]
          p2 = st.point_to_coords[point2]
          dp = p2-p1
          p = p1 + c1*dp
          rst.point_to_coords[rpoint] = p
          rp1 = st.point_to_rcoords[point1]
          rp2 = st.point_to_rcoords[point2]
          rdp = rp2-rp1
          rp = rp1 + c1*rdp
          rst.point_to_rcoords[rpoint] = rp
          for (ils, point_to_val) in enumerate(i_to_point_to_value)
            s1 = point_to_val[point1]
            s2 = point_to_val[point2]
            ds = s2-s1
            s = s1 + c1*ds
            i_to_rpoint_to_value[ils][rpoint] = s
          end
        end
      end
    end

    celloffset = rcell
    nsubcells = length(cell_table.case_to_subcell_to_inout[case])
    for subcell in 1:nsubcells
      rcell += 1
      rcell_to_inout[rcell] = cell_table.case_to_subcell_to_inout[case][subcell]
      rst.cell_to_bgcell[rcell] = st.cell_to_bgcell[cell]
      for (j,cell_to_inout) in enumerate(j_to_cell_to_inout)
        j_to_rcell_to_inout[j][rcell] = cell_to_inout[cell]
      end
      for subpoint in cell_table.case_to_subcell_to_points[case][subcell]
        q += 1
        rst.cell_to_points.data[q] = subpoint + pointoffset
      end
    end

  end
end

function initial_sub_triangulation(grid::Grid,geom::AnalyticalGeometry)
  geom_x = discretize(geom,grid)
  initial_sub_triangulation(grid,geom_x)
end

function initial_sub_triangulation(grid::Grid,geom::DiscreteGeometry)
  ugrid = UnstructuredGrid(grid)
  _initial_sub_triangulation(ugrid,geom)
end

function _initial_sub_triangulation(grid::UnstructuredGrid,geom::DiscreteGeometry)

  cutgrid, ls_to_cutpoint_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls = _extract_grid_of_cut_cells(grid,geom)

  subtrian, ls_to_subpoint_to_value = _simplexify_and_isolate_cells_in_cutgrid(cutgrid,ls_to_cutpoint_to_value)

  subtrian, ls_to_subpoint_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls
end

function _extract_grid_of_cut_cells(grid,geom)

  tree = get_tree(geom)
  ls_to_point_to_value, oid_to_ls = _find_unique_leaves(tree)

  p = _check_and_get_polytope(grid)
  table = LookupTable(p)
  cell_to_points = get_cell_nodes(grid)

  ls_to_cell_to_inoutcut = [
    _compute_in_out_or_cut(table,cell_to_points,point_to_value)
    for point_to_value in ls_to_point_to_value]

  cutcell_to_cell = _find_cut_cells(ls_to_cell_to_inoutcut)

  cutgrid = GridPortion(grid,cutcell_to_cell)

  ls_to_cutpoint_to_value = [
    point_to_value[cutgrid.node_to_oldnode] for point_to_value in ls_to_point_to_value ]

  cutgrid, ls_to_cutpoint_to_value, ls_to_cell_to_inoutcut, oid_to_ls
end

function _find_cut_cells(ls_to_cell_to_inoutcut)
  ncells = length(first(ls_to_cell_to_inoutcut))
  cell_to_iscut = fill(false,ncells)
  for cell in 1:ncells
    for cell_to_inoutcut in ls_to_cell_to_inoutcut
      inoutcut = cell_to_inoutcut[cell]
      if inoutcut == CUT
        cell_to_iscut[cell] = true
      end
    end
  end
  findall(cell_to_iscut)
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

function _simplexify_and_isolate_cells_in_cutgrid(cutgrid,ls_to_cutpoint_to_value)

  p = _check_and_get_polytope(cutgrid)

  ltcell_to_lpoints, simplex = simplexify(p)
  lpoint_to_lcoords = get_vertex_coordinates(p)
  _ensure_positive_jacobians!(ltcell_to_lpoints,lpoint_to_lcoords,simplex)

  out = _simplexify(
    get_node_coordinates(cutgrid),
    get_cell_nodes(cutgrid),
    ltcell_to_lpoints,
    lpoint_to_lcoords,
    ls_to_cutpoint_to_value,
    num_vertices(p),
    num_vertices(simplex))

  tcell_to_tpoints, tpoint_to_coords, tpoint_to_rcoords, ls_to_tpoint_to_value = out
  _ensure_positive_jacobians!(tcell_to_tpoints,tpoint_to_coords,simplex)

  ntcells = length(tcell_to_tpoints)
  nltcells = length(ltcell_to_lpoints)
  tcell_to_cell = _setup_cell_to_bgcell(cutgrid.cell_to_oldcell,nltcells,ntcells)

  subtrian = SubTriangulation(
    tcell_to_tpoints,
    tcell_to_cell,
    tpoint_to_coords,
    tpoint_to_rcoords)

  subtrian, ls_to_tpoint_to_value
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

