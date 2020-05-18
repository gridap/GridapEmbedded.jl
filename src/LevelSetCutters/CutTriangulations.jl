struct MiniCell
  num_points::Int
  edge_to_points::Vector{Vector{Int}}
end

function MiniCell(p::Polytope)
  MiniCell(num_vertices(p),get_faces(p,1,0))
end

function MiniCell(p::Polytope{0})
  MiniCell(1,[Int[]])
end

struct CutTriangulation{Dc,Dp,T}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dc,T}}
  cell_to_points::Table{Int,Int32}
  cell_to_bgcell::Vector{Int32}
  minicell::MiniCell
  table::LookupTable{Dc,T}
end

function CutTriangulation(
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dc,T}}
  cell_to_points::Table{Int,Int32}
  cell_to_bgcell::Vector{Int32},
  p::Polytope{Dc}) where {Dc,Dp,T}

  minicell = MiniCell(p)
  table = LookupTable(p)

  CutTriangulation(
    point_to_coords,
    point_to_rcoords,
    cell_to_points,
    cell_to_bgcell,
    minicell,
    table)
end

function allocate_sub_triangulation(
  m::CutTriangulation{Dc,Dp,T},
  n_cells::Integer,
  n_points::Integer) where {Dc,Dp,T}

  point_to_coords = zeros(Point{Dp,T},n_points)
  point_to_rcoords = zeros(Point{Dc,T},n_points)
  cell_to_points = _allocate_table(m.cell_to_points,n_cells,m.minicell.num_points)
  cell_to_bgcell = zeros(eltype(m.cell_to_bgcell),n_cells)
  cell_to_inoutcut = zeros(Int8,n_cells)

  s = CutTriangulation(
    point_to_coords,
    point_to_rcoords,
    cell_to_points,
    cell_to_bgcell,
    m.minicell,
    m.table)

  s, cell_to_inoutcut
end

function allocate_sub_triangulation_with_boundary(
  s::CutTriangulation{Dc,Dp,T},
  n_cells::Integer,
  n_points::Integer,
  n_facets::Integer) where {Dc,Dp,T}

  subcells, scell_to_inoutcut = allocate_sub_triangulation(s,n_cells,n_points)
  subfacets = allocate_boundary_triangulation(subcells,n_facets)
  subcells, scell_to_inoutcut, subfacets
end

function allocate_boundary_triangulation(
  s::CutTriangulation{Dc,Dp,T}, n_facets::Integer) where {Dc,Dp,T}

  f = Simplex(Val{Dc-1}())

  facet_to_points = _allocate_table(s.cell_to_points,n_facets,Dc)
  facet_to_bgcell = zeros(eltype(s.cell_to_bgcell),n_facets)
  facet_to_normal = zeros(VectorValue{Dp,T},n_facets)

  CutBoundaryTriangulation(
    s.point_to_coords,
    s.point_to_rcoords,
    facet_to_points,
    facet_to_bgcell,
    facet_to_normal,
    f)
end

function _allocate_table(a::Table{Td,Tp},n_cells,n_lpoints) where {Td,Tp}
  ndata = n_cells*n_lpoints
  data = zeros(Td,ndata)
  ptrs = fill!(Tp(n_lpoints),n_cells+1)
  length_to_ptrs!(ptrs)
  Table(data,ptrs)
end

get_minicell(a::CutTriangulation) = a.minicell

get_lookup_table(a::CutTriangulation) = a.table

get_cell_to_points(a::CutTriangulation) = a.cell_to_points

function set_cell_data!( s::CutTriangulation, m::CutTriangulation, d...)
  set_cell_data!(s.cell_to_bgcell, m.cell_to_bgcell, d...)
end

function set_point_data!( s::CutTriangulation, m::CutTriangulation, d...)
  set_point_data!(s.point_to_coords, m.point_to_coords, d...)
  set_point_data!(s.point_to_rcoords, m.point_to_rcoords, d...)
end

struct CutBoundaryTriangulation{Dc,Dp,T}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dc,T}}
  facet_to_points::Table{Int,Int32}
  facet_to_bgcell::Vector{Int32}
  facet_to_normal::Vector{Point{Dp,T}}
  minicell::MiniCell
  table::LookupTable{Dc,T}
end

function CutBoundaryTriangulation(
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dc,T}}
  facet_to_points::Table{Int,Int32}
  facet_to_bgcell::Vector{Int32},
  facet_to_normal::Vector{Point{Dp,T}},
  p::Polytope{Dc}) where {Dc,Dp,T}

  minicell = MiniCell(p)
  table = LookupTable(p)

  CutBoundaryTriangulation(
    point_to_coords,
    point_to_rcoords,
    facet_to_points,
    facet_to_bgcell,
    facet_to_normal,
    minicell,
    table)

end

function allocate_sub_triangulation(
  s::CutBoundaryTriangulation{Dc,Dp,T},
  n_facets::Integer,
  n_points::Integer) where {Dc,Dp,T}

  point_to_coords = zeros(Point{Dp,T},n_points)
  point_to_rcoords = zeros(Point{Dc,T},n_points)
  facet_to_points = _allocate_table(s.facet_to_points,n_facets,s.minicell.num_points)
  facet_to_bgcell = zeros(eltype(s.facet_to_bgcell),n_facets)
  facet_to_normal = zeros(VectorValue{Dp,T},n_facets)
  facet_to_inoutcut = zeros(Int8,n_facets)

  s = CutBoundaryTriangulation(
    point_to_coords,
    point_to_rcoords,
    facet_to_points,
    facet_to_bgcell,
    facet_to_normal,
    s.minicell,
    s.table)

  s, facet_to_inoutcut
end

get_minicell(a::CutBoundaryTriangulation) = a.minicell

get_lookup_table(a::CutBoundaryTriangulation) = a.table

get_cell_to_points(a::CutBoundaryTriangulation) = a.facet_to_points

get_point_to_coords(a::CutBoundaryTriangulation) = a.point_to_coords

function set_cell_data!( s::CutBoundaryTriangulation, m::CutBoundaryTriangulation, d...)
  set_cell_data!(s.facet_to_bgcell, m.facet_to_bgcell, d...)
  set_cell_data!(s.facet_to_normal, m.facet_to_normal, d...)
end

function set_cell_data!( s::CutBoundaryTriangulation, m::CutTriangulation, d...)
  set_cell_data!(s.facet_to_bgcell, m.cell_to_bgcell, d...)
end

function set_point_data!( s::CutBoundaryTriangulation, m::CutBoundaryTriangulation, d...)
  set_point_data!(s.point_to_coords, m.point_to_coords, d...)
  set_point_data!(s.point_to_rcoords, m.point_to_rcoords, d...)
end

struct CutTriangulatonWithLevelSet{A}
  sub_trian::A
  done_ls_to_cell_to_inoutcut::Vector{Vector{Int8}}
  pending_ls_to_point_to_value::Vector{Vector{Float64}}
end

function allocate_sub_triangulation(
  m::CutTriangulatonWithLevelSet,
  n_cells::Integer,
  n_points::Integer)

  s, scell_to_inoutcut = allocate_sub_triangulation(m.sub_trian,n_cells,n_points)

  n_done_ls = length(m.done_ls_to_cell_to_inoutcut)
  done_ls_to_cell_to_inoutcut = [ zeros(Int8,n_cells) for ls in 1:n_done_ls]

  n_pending_ls = length(m.pending_ls_to_point_to_value)
  pending_ls_to_point_to_value = [zeros(Float64,n_points) for ls in 1:n_pending_ls]

  _s = CutTriangulatonWithLevelSet(s,done_ls_to_cell_to_inoutcut,pending_ls_to_point_to_value)

  _s, scell_to_inoutcut
end

function allocate_sub_triangulation_with_boundary(
  s::CutTriangulatonWithLevelSet,
  n_cells::Integer,
  n_points::Integer,
  n_facets::Integer)

  s, scell_to_inoutcut = allocate_sub_triangulation(s,n_cells,n_points)
  s_boundary = allocate_boundary_triangulation(s.sub_trian,n_facets)
  s, scell_to_inoutcut, s_boundary
end

get_minicell(a::CutTriangulatonWithLevelSet) = get_minicell(a)

get_lookup_table(a::CutTriangulatonWithLevelSet) = get_lookup_table(a)

get_cell_to_points(a::CutTriangulatonWithLevelSet) = get_cell_to_points(a)

function set_cell_data!(s::CutTriangulatonWithLevelSet, m::CutTriangulatonWithLevelSet, d...)
  n_done_ls = length(s.done_ls_to_cell_to_inoutcut)
  for ls in  1:n_done_ls
    scell_to_inoutcut = s.done_ls_to_cell_to_inoutcut[ls]
    mcell_to_inoutcut = m.done_ls_to_cell_to_inoutcut[ls]
    set_cell_data!(scell_to_inoutcut, mcell_to_inoutcut, d...)
  end
  set_cell_data!(s.sub_trian, m.sub_trian, d...)
end

function set_cell_data!( s::CutBoundaryTriangulation, m::CutTriangulatonWithLevelSet, d...)
  set_cell_data!(s, m.sub_trian, d...)
end

function set_point_data!(s::CutTriangulatonWithLevelSet, m::CutTriangulatonWithLevelSet, d...)
  n_pending_ls = length(s.pending_ls_to_point_to_value)
  for ls in  1:n_pending_ls
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

function cut_sub_triangulation(m, mpoint_to_value)
  n_scells, n_spoints = count_sub_triangulation(m,point_to_value)
  s, scell_to_inoutcut  = allocate_sub_triangulation(m,n_scells,n_spoints)
  cut_sub_triangulation!(s,scell_to_inoutcut,m,mpoint_to_value)
  s, scell_to_inoutcut
end

function count_sub_triangulation(m,mpoint_to_value)
  n_scells = 0
  n_spoints = 0
  mcell_to_mpoints = get_cell_to_points(m)
  table = get_lookup_table(m)
  for mcell in 1:length(mcell_to_mpoints)
    case = compute_case(mcell_to_mpoints,mpoint_to_value,mcell)
    n_scells += length(table.case_to_subcell_to_inout[case])
    n_spoints += length(table.case_to_point_to_coordinates[case])
  end
  n_scells, n_spoints
end

function cut_sub_triangulation!(s,scell_to_inoutcut,m,mpoint_to_value)
  mcell_to_mpoints = get_cell_to_points(m)
  scell = 0
  spoint = 0
  q = 0
  for mcell in 1:length(mcell_to_mpoints)
    case = compute_case(mcell_to_mpoints,mpoint_to_value,mcell)
    scell, spoint, q = _cut_cell!(s,scell_to_inoutcut,m,mpoint_to_value,case,mcell,scell,spoint,q)
  end
end

function _cut_cell!(s,scell_to_inoutcut,m,mpoint_to_value,case,mcell,scell,spoint,q)

  mcell_to_mpoints = get_cell_to_points(m)
  scell_to_spoints = get_cell_to_points(s)
  minicell = get_minicell(m)
  table = get_lookup_table(m)

  nsubcells = length(table.case_to_subcell_to_inout[case])
  for subcell in 1:nsubcells
    scell += 1
    scell_to_inoutcut[scell] = table.case_to_subcell_to_inout[case][subcell]
    set_cell_data!(s,m,scell,mcell)
    for subpoint in table.case_to_subcell_to_points[case][subcell]
      q += 1
      scell_to_spoints.data[q] = subpoint + spoint
    end
  end

  a = mcell_to_mpoints.ptrs[mcell]-1

  for lpoint in 1:minicell.num_points
    spoint += 1
    mpoint = mcell_to_mpoints.data[a+lpoint]
    set_point_data!(s,m,spoint,mpoint)
  end

  if CUT == table.case_to_inoutcut[case]
    for (ledge, lpoints) in enumerate(minicell.edge_to_points)
      mpoint1 = mcell_to_mpoints.data[a+lpoints[1]]
      mpoint2 = mcell_to_mpoints.data[a+lpoints[2]]
      v1 = mpoint_to_value[mpoint1]
      v2 = mpoint_to_value[mpoint2]
      if isout(v1) != isout(v2)
        spoint += 1
        w1 = abs(v1)
        w2 = abs(v2)
        coeff = w1/(w1+w2)
        set_point_data!(s,m,spoint,mpoint1,mpoint2,coeff)
      end
    end
  end

  scell, spoint, q
end

function cut_sub_triangulation_with_boundary(m, mpoint_to_value)
  n_scells, n_spoints, n_facets = count_sub_triangulation_with_boundary(m,point_to_value)
  s, scell_to_inoutcut, s_boundary  = allocate_sub_triangulation_with_boundary(m,n_scells,n_spoints,n_facets)
  cut_sub_triangulation_with_boundary!(s,scell_to_inoutcut,s_boundary,m,mpoint_to_value)
  s, scell_to_inoutcut, s_boundary
end

function count_sub_triangulation_with_boundary(m,mpoint_to_value)
  n_scells = 0
  n_spoints = 0
  n_sfacets = 0
  mcell_to_mpoints = get_cell_to_points(m)
  table = get_lookup_table(m)
  for mcell in 1:length(mcell_to_mpoints)
    case = compute_case(mcell_to_mpoints,mpoint_to_value,mcell)
    n_scells += length(table.case_to_subcell_to_inout[case])
    n_spoints += length(table.case_to_point_to_coordinates[case])
    n_sfacets += length(table.case_to_subfacet_to_normal[case])
  end
  n_scells, n_spoints, n_sfacets
end

function cut_sub_triangulation_with_boundary!(s,scell_to_inoutcut,s_boundary,m,mpoint_to_value)
  mcell_to_mpoints = get_cell_to_points(m)
  scell = 0
  sfacet = 0
  spoint = 0
  q = 0
  z = 0
  for mcell in 1:length(mcell_to_mpoints)
    case = compute_case(mcell_to_mpoints,mpoint_to_value,mcell)
    pointoffset = spoint
    scell, spoint, q = _cut_cell!(s,scell_to_inoutcut,m,mpoint_to_value,case,mcell,scell,spoint,q)
    sfacet, z = _cut_cell_boundary!(s_boundary,m,case,mcell,sfacet,pointoffset,z)
  end
end

function  _cut_cell_boundary!(s_boundary,m,case,mcell,sfacet,spoint,z)

  table = get_lookup_table(m)
  sfacet_to_spoints = get_cell_to_points(s_boundary)
  sfacet_to_normal = s_boundary.facet_to_normal
  spoint_to_coords = get_point_to_coords(s_boundary)

  nsubfacets = length(table.case_to_subfacet_to_normal[case])
  for subfacet in 1:nsubfacets
    sfacet += 1
    for subpoint in table.case_to_subfacet_to_points[case][subfacet]
      z += 1
      sfacet_to_spoints.data[z] = subpoint + spoint
    end
    set_cell_data!(s_boundary,m,sfacet,mcell)
    normal = _setup_normal(
      table.case_to_subfacet_to_points[case],
      spoint_to_coords,
      subfacet,spoint)
    orientation = table.case_to_subfacet_to_orientation[case][subfacet]
    sfacet_to_normal[sfacet] = orientation*normal
  end

  sfacet, z
end

function cut_sub_triangulation_several_levelsets(sub_trian,pending_ls_to_point_to_value)

  pending_ls_to_point_to_value = copy(pending_ls_to_point_to_value)
  done_ls_to_cell_to_inoutcut = Vector{Int8}[]

  m = CutTriangulatonWithLevelSet(sub_trian,done_ls_to_cell_to_inoutcut,pending_ls_to_point_to_value)
  while length(m.pending_ls_to_point_to_value) > 0
    point_to_value = pop!(m.pending_ls_to_point_to_value)
    m, cell_to_inoutcut = cut_sub_triangulation(m, point_to_value)
    pushfirst!(m.done_ls_to_cell_to_inoutcut,cell_to_inoutcut)
  end

  m.sub_trian, m.done_ls_to_cell_to_inoutcut
end

function cut_sub_triangulation_with_boundary_several_levelsets(sub_trian,pending_ls_to_point_to_value)

  pending_ls_to_point_to_value = copy(pending_ls_to_point_to_value)
  done_ls_to_cell_to_inoutcut = Vector{Int8}[]

  m = CutTriangulatonWithLevelSet(sub_trian,done_ls_to_cell_to_inoutcut,pending_ls_to_point_to_value)
  while length(m.pending_ls_to_point_to_value) > 0


    point_to_value = pop!(m.pending_ls_to_point_to_value)
    m, cell_to_inoutcut, m_boundary = cut_sub_triangulation_with_boundary(m, point_to_value)

    m_boundary, done_ls_to_facet_to_inoutcut = cut_sub_triangulation_several_levelsets(
      m_boundary,m.pending_ls_to_point_to_value)


    pushfirst!(m.done_ls_to_cell_to_inoutcut,cell_to_inoutcut)
  end

  m.sub_trian, m.done_ls_to_cell_to_inoutcut

end

