## Initial sub-triangulation: IN/OUT/CUT of the background cells, simplices of the cut cells

function initial_sub_triangulation(grid::Grid,geom::AnalyticalGeometry)
  _, oid_to_ls = _find_unique_leaves(get_tree(geom))
  out = initial_sub_triangulation(grid,discretize(geom,grid))
  out[1], out[2], out[3], oid_to_ls
end

function initial_sub_triangulation(grid::Grid,geom::DiscreteGeometry)
  ugrid = UnstructuredGrid(grid)
  ls_to_point_to_value, oid_to_ls = _find_unique_leaves(get_tree(geom))
  out = _initial_sub_triangulation(ugrid,ls_to_point_to_value)
  out[1], out[2], out[3], oid_to_ls
end

function _initial_sub_triangulation(grid::UnstructuredGrid,ls_to_point_to_value)
  cutgrid, ls_to_cutpoint_to_value, ls_to_bgcell_to_inoutcut =
    _extract_grid_of_cut_cells(grid,ls_to_point_to_value)
  subtrian, ls_to_subpoint_to_value =
    _simplexify_and_isolate_cells_in_cutgrid(cutgrid,ls_to_cutpoint_to_value)
  subtrian, ls_to_subpoint_to_value, ls_to_bgcell_to_inoutcut
end

function _extract_grid_of_cut_cells(grid,ls_to_point_to_value)
  cell_to_points = get_cell_node_ids(grid)

  ls_to_cell_to_inoutcut = [
    _compute_in_out_or_cut(cell_to_points,point_to_value)
    for point_to_value in ls_to_point_to_value]

  cutcell_to_cell = _find_cut_cells(ls_to_cell_to_inoutcut)

  cutgrid = GridPortion(grid,cutcell_to_cell)

  ls_to_cutpoint_to_value = [
    point_to_value[cutgrid.node_to_parent_node] for point_to_value in ls_to_point_to_value ]

  return cutgrid, ls_to_cutpoint_to_value, ls_to_cell_to_inoutcut
end

# IN/OUT/CUT of every cell from the values of the level set at its nodes (no lookup table
# needed: the rule of `_find_in_out_or_cut`, the same one the tables are built with)
function _compute_in_out_or_cut(cell_to_points::Table,point_to_value::AbstractVector)
  ncells = length(cell_to_points)
  cell_to_inoutcut = zeros(Int8,ncells)
  for cell in 1:ncells
    cell_to_inoutcut[cell] = _find_in_out_or_cut(view(point_to_value,dataview(cell_to_points,cell)))
  end
  return cell_to_inoutcut
end

# Cells cut by at least one level set
function _find_cut_cells(ls_to_cell_to_inoutcut)
  ncells = length(first(ls_to_cell_to_inoutcut))
  findall(cell -> any(c -> c[cell] == CUT, ls_to_cell_to_inoutcut), 1:ncells)
end

# The polytope of a grid with a single type of first-order cells
function _check_and_get_polytope(grid)
  reffes = get_reffes(grid)
  @notimplementedif length(reffes) != 1
  reffe = first(reffes)
  @notimplementedif get_order(reffe) != 1
  return get_polytope(reffe)
end

# Splits every cut cell into simplices with their own (duplicated) points
function _simplexify_and_isolate_cells_in_cutgrid(cutgrid,ls_to_cutpoint_to_value)

  Dc = num_cell_dims(cutgrid)
  Dp = num_point_dims(cutgrid)

  p = _check_and_get_polytope(cutgrid)

  ltcell_to_lpoints, simplex = simplexify(p)
  lpoint_to_lcoords = get_vertex_coordinates(p)
  _ensure_positive_jacobians!(ltcell_to_lpoints,lpoint_to_lcoords,simplex)

  out = _simplexify(
    get_node_coordinates(cutgrid),
    get_cell_node_ids(cutgrid),
    ltcell_to_lpoints,
    lpoint_to_lcoords,
    ls_to_cutpoint_to_value,
    num_vertices(p),
    num_vertices(simplex),
    Dc)

  tcell_to_tpoints, tpoint_to_coords, tpoint_to_rcoords, ls_to_tpoint_to_value = out

  nltcells = length(ltcell_to_lpoints)
  tcell_to_cell = _setup_cell_to_bgcell(cutgrid.cell_to_parent_cell,nltcells)

  if Dc == Dp
    _ensure_positive_jacobians!(tcell_to_tpoints,tpoint_to_coords,simplex)
  elseif Dc+1 == Dp
    tcell_to_cutcell = _setup_cell_to_bgcell(1:num_cells(cutgrid),nltcells)
    _ensure_positive_jacobians_facets!(tcell_to_tpoints,tpoint_to_coords,simplex,cutgrid,tcell_to_cutcell)
  else
    @notimplemented
  end

  subtrian = SubCellData(
    tcell_to_tpoints,
    tcell_to_cell,
    tpoint_to_coords,
    tpoint_to_rcoords)

  return subtrian, ls_to_tpoint_to_value
end

# Each parent cell repeated for its `nlcells` simplices
function _setup_cell_to_bgcell(pcell_to_bgcell,nlcells)
  Int32[ bgcell for bgcell in pcell_to_bgcell for lcell in 1:nlcells ]
end

function _simplexify(
  point_to_coords,
  cell_to_points::Table,
  ltcell_to_lpoints,
  lpoint_to_lcoords,
  ls_to_point_to_value,
  nlpoints,
  nsp,
  Dc)

  ncells = length(cell_to_points)
  nltcells = length(ltcell_to_lpoints)
  ntcells = ncells*nltcells
  ntpoints = ncells*nlpoints
  T = eltype(first(ls_to_point_to_value))
  tcell_to_tpoints_data = zeros(eltype(cell_to_points.data),nsp*ntcells)
  tcell_to_tpoints_ptrs = fill(eltype(cell_to_points.ptrs)(nsp),ntcells+1)
  length_to_ptrs!(tcell_to_tpoints_ptrs)
  tcell_to_tpoints = Table(tcell_to_tpoints_data,tcell_to_tpoints_ptrs)
  Dp = eltype(point_to_coords).parameters[1]
  tpoint_to_coords = zeros(Point{Dp,T},ntpoints)
  tpoint_to_rcoords = zeros(Point{Dc,T},ntpoints)
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

  return tcell_to_tpoints, tpoint_to_coords, tpoint_to_rcoords, ls_to_tpoint_to_value
end

function _ensure_positive_jacobians_facets!(
  tcell_to_tpoints,tpoint_to_coords,simplex,cutgrid,tcell_to_cutcell)

  n_tcells = length(tcell_to_tpoints)
  tcell_to_ctype = ones(Int8,n_tcells)
  ctype_to_reffe = [LagrangianRefFE(Float64,simplex,1)]
  tgrid = UnstructuredGrid(tpoint_to_coords,tcell_to_tpoints,ctype_to_reffe,tcell_to_ctype)

  tmap = get_cell_map(tgrid)
  cutmap = lazy_map(Reindex(get_cell_map(cutgrid)),tcell_to_cutcell)

  ctype_to_q = map(r->[first(get_node_coordinates(r))],ctype_to_reffe)
  tcell_to_q = expand_cell_data(ctype_to_q,tcell_to_ctype)

  tjac = lazy_map(∇,tmap)
  jac = lazy_map(∇,cutmap)
  tjac_q = lazy_map(evaluate,tjac,tcell_to_q)
  jac_q = lazy_map(evaluate,jac,tcell_to_q)

  c1 = array_cache(tjac_q)
  c2 = array_cache(jac_q)

  _ensure_positive_jacobians_facets_work!(tcell_to_tpoints,c1,c2,tjac_q,jac_q)
  return nothing
end

function _ensure_positive_jacobians_facets_work!(tcell_to_tpoints,c1,c2,tjac_q,jac_q)
  n_tcells = length(tcell_to_tpoints)
  for tcell in 1:n_tcells
    tjac = getindex!(c1,tjac_q,tcell)
    jac = getindex!(c2,jac_q,tcell)
    n1 = VectorValue(first(tjac)...)
    n2 = VectorValue(first(jac)...)
    if n1⋅n2 < 0
      points = dataview(tcell_to_tpoints,tcell)
      points[1], points[2] = points[2], points[1]
    end
  end
  return nothing
end

## Cutting by several level sets, one after the other

function cut_sub_triangulation_several_levelsets(sub_trian,pending_ls_to_point_to_value)
  pending_ls_to_point_to_value = copy(pending_ls_to_point_to_value)
  done_ls_to_cell_to_inoutcut = Vector{Int8}[]

  m = CutTriangulation(sub_trian,done_ls_to_cell_to_inoutcut,pending_ls_to_point_to_value)
  while length(m.pending_ls_to_point_to_value) > 0
    point_to_value = pop!(m.pending_ls_to_point_to_value)
    m, cell_to_inoutcut = cut_sub_triangulation(m, point_to_value)
    pushfirst!(m.done_ls_to_cell_to_inoutcut,cell_to_inoutcut)
  end

  return m.sub_trian, m.done_ls_to_cell_to_inoutcut
end

function cut_sub_triangulation_with_boundary_several_levelsets(sub_trian,pending_ls_to_point_to_value)
  pending_ls_to_point_to_value = copy(pending_ls_to_point_to_value)
  done_ls_to_cell_to_inoutcut = Vector{Int8}[]
  done_ls_to_boundary = []
  done_ls_to_ls_to_facet_inoutcut = Vector{Vector{Int8}}[]

  m = CutTriangulation(sub_trian,done_ls_to_cell_to_inoutcut,pending_ls_to_point_to_value)

  n_ls = length(pending_ls_to_point_to_value)
  for ls in n_ls:-1:1

    point_to_value = m.pending_ls_to_point_to_value[ls]
    m, cell_to_inoutcut, m_boundary = cut_sub_triangulation_with_boundary(m, point_to_value)

    ls_to_point_to_value = deleteat!(copy(m.pending_ls_to_point_to_value),ls)

    m_boundary, ls_to_facet_to_inoutcut =
      cut_sub_triangulation_several_levelsets(m_boundary,ls_to_point_to_value)

    n_facets = length(get_cell_to_points(m_boundary))
    facet_to_inoutcut = fill(Int8(INTERFACE),n_facets)
    insert!(ls_to_facet_to_inoutcut,ls,facet_to_inoutcut)

    pushfirst!(m.done_ls_to_cell_to_inoutcut,cell_to_inoutcut)
    pushfirst!(done_ls_to_boundary,m_boundary)
    pushfirst!(done_ls_to_ls_to_facet_inoutcut,ls_to_facet_to_inoutcut)
  end

  m_boundary, done_ls_to_facet_to_inoutcut =
    merge_sub_face_data(done_ls_to_boundary,done_ls_to_ls_to_facet_inoutcut)

  return m.sub_trian, m.done_ls_to_cell_to_inoutcut, m_boundary, done_ls_to_facet_to_inoutcut
end

## A sub-triangulation being cut: level sets already applied (labels) and pending (values)

struct CutTriangulation{Dc,A,T}
  sub_trian::A
  done_ls_to_cell_to_inoutcut::Vector{Vector{Int8}}
  pending_ls_to_point_to_value::Vector{Vector{T}}
  table::LookupTable{Dc}
end

function CutTriangulation(
  sub_trian,
  done_ls_to_cell_to_inoutcut::Vector{Vector{Int8}},
  pending_ls_to_point_to_value::Vector{Vector{T}}) where T

  Dc = get_cell_dim(sub_trian)
  table = get_lookup_table(Simplex(Val{Dc}()))
  return CutTriangulation(
    sub_trian, done_ls_to_cell_to_inoutcut, pending_ls_to_point_to_value, table)
end

get_lookup_table(a::CutTriangulation) = a.table

get_cell_to_points(a::CutTriangulation) = get_cell_to_points(a.sub_trian)

get_facet_to_normal(a::CutTriangulation) = get_facet_to_normal(a.sub_trian)

get_point_to_coords(a::CutTriangulation) = get_point_to_coords(a.sub_trian)

function allocate_sub_triangulation(
  m::CutTriangulation, n_cells::Integer, n_points::Integer)

  sub_trian, scell_to_inoutcut = allocate_sub_triangulation(m.sub_trian,n_cells,n_points)

  n_done_ls = length(m.done_ls_to_cell_to_inoutcut)
  done_ls_to_cell_to_inoutcut = [ zeros(Int8,n_cells) for ls in 1:n_done_ls]

  n_pending_ls = length(m.pending_ls_to_point_to_value)
  pending_ls_to_point_to_value = [zeros(eltype(eltype(m.pending_ls_to_point_to_value)),n_points) for ls in 1:n_pending_ls]

  s = CutTriangulation(
    sub_trian,done_ls_to_cell_to_inoutcut,pending_ls_to_point_to_value)

  return s, scell_to_inoutcut
end

function allocate_sub_triangulation_with_boundary(
  m::CutTriangulation,
  n_cells::Integer,
  n_points::Integer,
  n_facets::Integer)

  s, scell_to_inoutcut = allocate_sub_triangulation(m,n_cells,n_points)
  sub_trian_boundary = allocate_boundary_triangulation(s.sub_trian,n_facets)

  return s, scell_to_inoutcut, sub_trian_boundary
end

function set_cell_data!(s::CutTriangulation, m::CutTriangulation, d...)
  n_done_ls = length(s.done_ls_to_cell_to_inoutcut)
  for ls in 1:n_done_ls
    scell_to_inoutcut = s.done_ls_to_cell_to_inoutcut[ls]
    mcell_to_inoutcut = m.done_ls_to_cell_to_inoutcut[ls]
    set_cell_data!(scell_to_inoutcut, mcell_to_inoutcut, d...)
  end
  set_cell_data!(s.sub_trian, m.sub_trian, d...)
end

function set_cell_data_boundary!(s, m::CutTriangulation, d...)
  set_cell_data_boundary!(s, m.sub_trian, d...)
end

function set_point_data!(s::CutTriangulation, m::CutTriangulation, d...)
  n_pending_ls = length(s.pending_ls_to_point_to_value)
  for ls in 1:n_pending_ls
    spoint_to_value = s.pending_ls_to_point_to_value[ls]
    mpoint_to_value = m.pending_ls_to_point_to_value[ls]
    set_point_data!(spoint_to_value, mpoint_to_value, d...)
  end
  set_point_data!(s.sub_trian, m.sub_trian, d...)
end

function set_cell_data!(
  scell_to_value::AbstractVector,
  mcell_to_value::AbstractVector,
  scell::Integer,
  mcell::Integer)

  mvalue = mcell_to_value[mcell]
  scell_to_value[scell] = mvalue
end

function set_point_data!(
  spoint_to_value::AbstractVector,
  mpoint_to_value::AbstractVector,
  spoint::Integer,
  mpoint::Integer)

  mvalue = mpoint_to_value[mpoint]
  spoint_to_value[spoint] = mvalue
end

function set_point_data!(
  spoint_to_value::AbstractVector,
  mpoint_to_value::AbstractVector,
  spoint::Integer,
  mpoint1::Integer,
  mpoint2::Integer,
  coeff::Real)

  mvalue1 = mpoint_to_value[mpoint1]
  mvalue2 = mpoint_to_value[mpoint2]
  svalue = mvalue1 + coeff*(mvalue2 - mvalue1)
  spoint_to_value[spoint] = svalue
end

## Cutting by one level set (count, allocate, fill)

function cut_sub_triangulation(m, mpoint_to_value)
  _m = CutTriangulation(m,Vector{Int8}[],Vector{Float64}[])
  _s, cell_to_inoutcut = cut_sub_triangulation(_m,mpoint_to_value)
  s = _s.sub_trian
  return s, cell_to_inoutcut
end

function cut_sub_triangulation(m::CutTriangulation, mpoint_to_value)
  n_scells, n_spoints, _ = count_sub_triangulation(m,mpoint_to_value)
  s, scell_to_inoutcut  = allocate_sub_triangulation(m,n_scells,n_spoints)
  cut_sub_triangulation!(s,scell_to_inoutcut,m,mpoint_to_value)
  return s, scell_to_inoutcut
end

function cut_sub_triangulation!(s,scell_to_inoutcut,m,mpoint_to_value)
  mcell_to_mpoints = get_cell_to_points(m)
  scell, spoint, q = 0, 0, 0
  for mcell in 1:length(mcell_to_mpoints)
    case = compute_case(mcell_to_mpoints,mpoint_to_value,mcell)
    scell, spoint, q = _cut_cell!(s,scell_to_inoutcut,m,mpoint_to_value,case,mcell,scell,spoint,q)
  end
end

function cut_sub_triangulation_with_boundary(m, mpoint_to_value)
  _m = CutTriangulation(m,Vector{Int8}[],Vector{Float64}[])
  _s, cell_to_inoutcut, s_boundary = cut_sub_triangulation_with_boundary(_m,mpoint_to_value)
  s = _s.sub_trian
  return s, cell_to_inoutcut, s_boundary
end

function cut_sub_triangulation_with_boundary(m::CutTriangulation, mpoint_to_value)
  n_scells, n_spoints, n_facets = count_sub_triangulation(m,mpoint_to_value)
  s, scell_to_inoutcut, s_boundary  = allocate_sub_triangulation_with_boundary(m,n_scells,n_spoints,n_facets)
  cut_sub_triangulation_with_boundary!(s,scell_to_inoutcut,s_boundary,m,mpoint_to_value)
  return s, scell_to_inoutcut, s_boundary
end

function cut_sub_triangulation_with_boundary!(s,scell_to_inoutcut,s_boundary,m,mpoint_to_value)
  q, z = 0, 0
  scell, sfacet, spoint = 0, 0, 0
  mcell_to_mpoints = get_cell_to_points(m)
  for mcell in eachindex(mcell_to_mpoints)
    case = compute_case(mcell_to_mpoints,mpoint_to_value,mcell)
    pointoffset = spoint
    scell, spoint, q = _cut_cell!(s,scell_to_inoutcut,m,mpoint_to_value,case,mcell,scell,spoint,q)
    sfacet, z = _cut_cell_boundary!(s_boundary,m,case,mcell,sfacet,pointoffset,z)
  end
end

# Number of subcells, points and subfacets produced by cutting `m` with `mpoint_to_value`
function count_sub_triangulation(m,mpoint_to_value)
  mcell_to_mpoints = get_cell_to_points(m)
  table = get_lookup_table(m)
  n_scells = 0
  n_spoints = 0
  n_sfacets = 0
  for mcell in 1:length(mcell_to_mpoints)
    case = compute_case(mcell_to_mpoints,mpoint_to_value,mcell)
    n_scells += length(table.case_to_subcell_to_points[case])
    n_spoints += num_case_points(table,case)
    n_sfacets += length(table.case_to_subfacet_to_points[case])
  end
  n_scells, n_spoints, n_sfacets
end

# Writes the subcells and points of `mcell` after the `scell` cells, `spoint` points and `q`
# connectivity entries already written; returns the updated counters.
function _cut_cell!(s,scell_to_inoutcut,m,mpoint_to_value,case,mcell,scell,spoint,q)

  mcell_to_mpoints = get_cell_to_points(m)
  scell_to_spoints = get_cell_to_points(s)
  table = get_lookup_table(m)

  # Subcells: the table's connectivity, shifted by the points already written
  subcell_to_inout = table.case_to_subcell_to_inout[case]
  for (subcell,subpoints) in enumerate(table.case_to_subcell_to_points[case])
    scell += 1
    scell_to_inoutcut[scell] = subcell_to_inout[subcell]
    set_cell_data!(s,m,scell,mcell)
    for subpoint in subpoints
      q += 1
      scell_to_spoints.data[q] = subpoint + spoint
    end
  end

  # Points, in the table's order: the vertices of the cell, then one per cut edge
  # (linear interpolation of the zero of the level set along the edge)
  a = mcell_to_mpoints.ptrs[mcell]-1
  for lvertex in 1:table.num_vertices
    spoint += 1
    mpoint = mcell_to_mpoints.data[a+lvertex]
    set_point_data!(s,m,spoint,mpoint)
  end
  for (lvertex1,lvertex2) in table.case_to_cut_edges[case]
    mpoint1 = mcell_to_mpoints.data[a+lvertex1]
    mpoint2 = mcell_to_mpoints.data[a+lvertex2]
    w1 = abs(mpoint_to_value[mpoint1])
    w2 = abs(mpoint_to_value[mpoint2])
    coeff = w1/(w1+w2)
    spoint += 1
    set_point_data!(s,m,spoint,mpoint1,mpoint2,coeff)
  end

  return scell, spoint, q
end

# Writes the interface pieces of `mcell`; `spoint` is the offset of the cell's points
function _cut_cell_boundary!(s_boundary,m,case,mcell,sfacet,spoint,z)
  table = get_lookup_table(m)
  sfacet_to_spoints = get_cell_to_points(s_boundary)
  sfacet_to_normal = get_facet_to_normal(s_boundary)
  spoint_to_coords = get_point_to_coords(s_boundary)

  subfacet_to_points = table.case_to_subfacet_to_points[case]
  for (subfacet,subpoints) in enumerate(subfacet_to_points)
    sfacet += 1
    for subpoint in subpoints
      z += 1
      sfacet_to_spoints.data[z] = subpoint + spoint
    end
    set_cell_data_boundary!(s_boundary,m,sfacet,mcell)
    normal = _setup_normal(subfacet_to_points,spoint_to_coords,subfacet,spoint)
    orientation = table.case_to_subfacet_to_orientation[case][subfacet]
    sfacet_to_normal[sfacet] = orientation*normal
  end

  return sfacet, z
end

## SubFaceData (SubCellData / SubFacetData) as the sub-triangulation being cut

function allocate_sub_triangulation(
  m::SubCellData{Dc,Dp,T}, n_cells::Integer, n_points::Integer) where {Dc,Dp,T}

  point_to_coords = zeros(Point{Dp,T},n_points)
  point_to_rcoords = zeros(Point{Dc,T},n_points)
  cell_to_points = _allocate_table(m.cell_to_points,n_cells,Dc+1)
  cell_to_bgcell = zeros(eltype(m.cell_to_bgcell),n_cells)
  cell_to_inoutcut = zeros(Int8,n_cells)

  s = SubCellData(
    cell_to_points, cell_to_bgcell, point_to_coords, point_to_rcoords)

  s, cell_to_inoutcut
end

function allocate_boundary_triangulation(
  s::SubCellData{Dc,Dp,T}, n_facets::Integer) where {Dc,Dp,T}

  facet_to_points = _allocate_table(s.cell_to_points,n_facets,Dc)
  facet_to_bgcell = zeros(eltype(s.cell_to_bgcell),n_facets)
  facet_to_normal = zeros(VectorValue{Dp,T},n_facets)

  SubFacetData(
    facet_to_points,
    facet_to_normal,
    facet_to_bgcell,
    s.point_to_coords,
    s.point_to_rcoords)
end

function _allocate_table(a::Table{T,Vd,Vp},n_cells,n_lpoints) where {T,Vd,Vp}
  ndata = n_cells*n_lpoints
  data = zeros(T,ndata)
  ptrs = fill(eltype(Vp)(n_lpoints),n_cells+1)
  length_to_ptrs!(ptrs)
  Table(data,ptrs)
end

get_cell_to_points(m::SubFaceData) = m.cell_to_points

get_cell_dim(m::SubFaceData{Df}) where Df = Df

get_point_to_coords(m::SubFaceData) = m.point_to_coords

function set_point_data!(s::SubFaceData, m::SubFaceData, d...)
  set_point_data!(s.point_to_coords, m.point_to_coords, d...)
  set_point_data!(s.point_to_rcoords, m.point_to_rcoords, d...)
end

function set_cell_data!(s::SubCellData, m::SubCellData, d...)
  set_cell_data!(s.cell_to_bgcell, m.cell_to_bgcell, d...)
end

function set_cell_data_boundary!(s::SubFacetData, m::SubCellData, d...)
  set_cell_data!(s.cell_to_bgcell, m.cell_to_bgcell, d...)
end

function allocate_sub_triangulation(
  m::SubFacetData{Dp,T}, n_facets::Integer, n_points::Integer) where {Dp,T}

  point_to_coords = zeros(Point{Dp,T},n_points)
  point_to_rcoords = zeros(Point{Dp,T},n_points)
  facet_to_points = _allocate_table(m.cell_to_points,n_facets,Dp)
  facet_to_bgcell = zeros(eltype(m.cell_to_bgcell),n_facets)
  facet_to_normal = zeros(VectorValue{Dp,T},n_facets)
  facet_to_inoutcut = zeros(Int8,n_facets)

  s = SubFacetData(
    facet_to_points,
    facet_to_normal,
    facet_to_bgcell,
    point_to_coords,
    point_to_rcoords)

  s, facet_to_inoutcut
end

get_facet_to_normal(m::SubFacetData) = m.cell_to_normal

function set_cell_data!(s::SubFacetData, m::SubFacetData, d...)
  set_cell_data!(s.cell_to_bgcell, m.cell_to_bgcell, d...)
  set_cell_data!(s.cell_to_normal, m.cell_to_normal, d...)
end
