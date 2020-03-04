
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

function initial_sub_triangulation(_grid::Grid,_ls_to_point_to_value::Vector{<:AbstractArray})
  _ls_to_point_to_value_ = map(collect1d,_ls_to_point_to_value)
  grid = UnstructuredGrid(_grid)
  reffes = get_reffes(grid)
  @notimplementedif length(reffes) != 1
  reffe = first(reffes)
  order = 1
  @notimplementedif get_order(reffe) != order
  p = get_polytope(reffe)
  _table = LookupTable(p)
  _cell_to_points = get_cell_nodes(grid)
  _cell_to_inoutcut = compute_in_out_or_cut(_table,_cell_to_points,_ls_to_point_to_value_)
  cutcell_to_cell = findall(_cell_to_inoutcut .== CUT)
  cutgrid = GridPortion(grid,cutcell_to_cell)
  ls_to_point_to_value = [ point_to_value[cutgrid.node_to_oldnode] for point_to_value in _ls_to_point_to_value_ ]
  ltcell_to_lpoints, simplex = simplexify(p)
  lpoint_to_lcoords = get_vertex_coordinates(p)
  nlpoints = num_vertices(p)
  nsp = num_vertices(simplex)

  tcell_to_tpoints, tpoint_to_coords, tpoint_to_rcoords, ls_to_tpoint_to_value = _simplexity(
    get_node_coordinates(cutgrid),
    get_cell_nodes(cutgrid),
    ltcell_to_lpoints,
    lpoint_to_lcoords,
    ls_to_point_to_value,
    nlpoints,
    nsp)

  ntcells = length(tcell_to_tpoints)
  tcell_to_inoutcut = fill(Int8(CUT),ntcells)
  nltcells = length(ltcell_to_lpoints)
  tcell_to_bgcell = _setup_cell_to_bgcell(cutgrid.cell_to_oldcell,nltcells,ntcells)
  table = LookupTable(simplex)

  SubTriangulation(
    table,
    tcell_to_tpoints,
    tcell_to_inoutcut,
    tcell_to_bgcell,
    tpoint_to_coords,
    ls_to_tpoint_to_value)
end

function _simplexity(
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

function cut_sub_triangulation(st::SubTriangulation)
  _st = st
  while length(_st.ls_to_point_to_value) > 0
    point_to_value = pop!(_st.ls_to_point_to_value)
    _st = _cut_sub_triangulation(_st,point_to_value)
  end
  _st
end

function _cut_sub_triangulation(st::SubTriangulation,point_to_value)
  nrcells, nrpoints = _cut_sub_triangulation_count(st,point_to_value)
  rst = _allocate_new_sub_triangulation(st,nrcells,nrpoints)
  _fill_new_sub_triangulation!(rst,st,point_to_value)
  rst
end

function _cut_sub_triangulation_count(st,point_to_value)
  nrcells = 0
  nrpoints = 0
  for cell in 1:length(st.cell_to_points)
    case = compute_case(st.cell_to_points,point_to_value,cell)
    nrcells += length(st.table.case_to_subcell_to_inout[case])
    nrpoints += length(st.table.case_to_point_to_coordinates[case])
  end
  nrcells, nrpoints
end

function _allocate_new_sub_triangulation(st::SubTriangulation{D,T},nrcells,nrpoints) where {D,T}
  nlp = st.table.nlpoints
  rcell_to_rpoints_data = zeros(eltype(st.cell_to_points.data),nlp*nrcells)
  rcell_to_rpoints_ptrs = fill(eltype(st.cell_to_points.ptrs)(nlp),nrcells+1)
  length_to_ptrs!(rcell_to_rpoints_ptrs)
  rcell_to_rpoints = Table(rcell_to_rpoints_data,rcell_to_rpoints_ptrs)
  rcell_to_inoutcut = zeros(eltype(st.cell_to_inoutcut),nrcells)
  rcell_to_bgcell = zeros(eltype(st.cell_to_bgcell),nrcells)
  rpoint_to_coords = zeros(Point{D,T},nrpoints)
  ls_to_rpoint_to_value = [zeros(T,nrpoints) for i in 1:length(st.ls_to_point_to_value)]

  SubTriangulation(
    st.table,
    rcell_to_rpoints,
    rcell_to_inoutcut,
    rcell_to_bgcell,
    rpoint_to_coords,
    ls_to_rpoint_to_value)

end

function _fill_new_sub_triangulation!(rst,st,point_to_value)
  rcell = 0
  rpoint = 0
  q = 0
  for cell in 1:length(st.cell_to_points)

    offset = rpoint
    a = st.cell_to_points.ptrs[cell]-1
    for lpoint in 1:st.table.nlpoints
      rpoint += 1
      point = st.cell_to_points.data[a+lpoint]
      coords = st.point_to_coords[point]
      rst.point_to_coords[rpoint] = coords
      for (ils, point_to_val) in enumerate(st.ls_to_point_to_value)
        val = point_to_val[point]
        rst.ls_to_point_to_value[ils][rpoint] = val
      end
    end

    case = compute_case(st.cell_to_points,point_to_value,cell)

    if CUT == st.table.case_to_inoutcut[case]
      for (ledge, lpoints) in enumerate(st.table.ledge_to_lpoints)
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
          for (ils, point_to_val) in enumerate(st.ls_to_point_to_value)
            s1 = point_to_val[point1]
            s2 = point_to_val[point2]
            ds = s2-s1
            s = s1 + c1*ds
            rst.ls_to_point_to_value[ils][rpoint] = s
          end
        end
      end
    end

    nsubcells = length(st.table.case_to_subcell_to_inout[case])
    for subcell in 1:nsubcells
      rcell += 1
      if st.cell_to_inoutcut[cell] == OUT
        rst.cell_to_inoutcut[rcell] = OUT
      else
        rst.cell_to_inoutcut[rcell] = st.table.case_to_subcell_to_inout[case][subcell]
      end
      rst.cell_to_bgcell[rcell] = st.cell_to_bgcell[cell]
      for subpoint in st.table.case_to_subcell_to_points[case][subcell]
        q += 1
        rst.cell_to_points.data[q] = subpoint + offset
      end
    end

  end
end

