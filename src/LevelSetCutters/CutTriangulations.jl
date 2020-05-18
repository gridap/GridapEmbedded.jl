
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

  CutTriangulation(
    point_to_coords,
    point_to_rcoords,
    cell_to_points,
    cell_to_bgcell,
    m.minicell,
    m.table)
end

function allocate_sub_triangulation_with_boundary(
  s::CutTriangulation{Dc,Dp,T},
  n_cells::Integer,
  n_points::Integer,
  n_facets::Integer) where {Dc,Dp,T}

  subcells = allocate_sub_triangulation(s,n_cells,n_points)

  f = Simplex(Val{Dc-1}())

  facet_to_points = _allocate_table(s.cell_to_points,n_facets,Dc)
  facet_to_bgcell = zeros(eltype(s.cell_to_bgcell),n_facets)
  facet_to_normal = zeros(VectorValue{Dp,T},n_facets)

  subfacets = CutBoundaryTriangulation(
    subcells.point_to_coords,
    subcells.point_to_rcoords,
    facet_to_points,
    facet_to_bgcell,
    facet_to_normal,
    f)

  subcells, subfacets

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

  CutBoundaryTriangulation(
    point_to_coords,
    point_to_rcoords,
    facet_to_points,
    facet_to_bgcell,
    facet_to_normal,
    s.minicell,
    s.table)
end

get_minicell(a::CutBoundaryTriangulation) = a.minicell

get_lookup_table(a::CutBoundaryTriangulation) = a.table

get_cell_to_points(a::CutBoundaryTriangulation) = a.facet_to_points

function set_cell_data!( s::CutBoundaryTriangulation, m::CutBoundaryTriangulation, d...)
  set_cell_data!(s.facet_to_bgcell, m.facet_to_bgcell, d...)
  set_cell_data!(s.facet_to_normal, m.facet_to_normal, d...)
end

function set_point_data!( s::CutBoundaryTriangulation, m::CutBoundaryTriangulation, d...)
  set_point_data!(s.point_to_coords, m.point_to_coords, d...)
  set_point_data!(s.point_to_rcoords, m.point_to_rcoords, d...)
end

struct SubtriangulatonWithLevelSet{A}
  sub_trian::A
  done_ls_to_cell_to_inoutcut::Vector{Vector{Int8}}
  pending_ls_to_point_to_value::Vector{Vector{Int8}}
end

get_minicell(a::SubtriangulatonWithLevelSet) = get_minicell(a)

get_lookup_table(a::SubtriangulatonWithLevelSet) = get_lookup_table(a)

get_cell_to_points(a::SubtriangulatonWithLevelSet) = get_cell_to_points(a)

function set_cell_data!(s::SubtriangulatonWithLevelSet, m::SubtriangulatonWithLevelSet, d...)
  n_done_ls = length(s.done_ls_to_cell_to_inoutcut)
  for ls in  1:n_done_ls
    scell_to_inoutcut = s.done_ls_to_cell_to_inoutcut[ls]
    mcell_to_inoutcut = m.done_ls_to_cell_to_inoutcut[ls]
    set_cell_data!(scell_to_inoutcut, mcell_to_inoutcut, d...)
  end
  set_cell_data!(s.sub_trian, m.sub_trian, d...)
end

function set_point_data!(s::SubtriangulatonWithLevelSet, m::SubtriangulatonWithLevelSet, d...)
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

function cut_sub_triangulation(m::SubtriangulatonWithLevelSet, mpoint_to_value::AbstractVector)

  n_scells, n_spoints = _cut_sub_triangulation_count(m,point_to_value)
  s, scell_to_inoutcut  = allocate_sub_triangulation(m,n_scells,n_spoints)

  s, scell_to_inoutcut
end

function cut_sub_triangulation_with_boundary(m::SubtriangulatonWithLevelSet, mpoint_to_value::AbstractVector)

  n_scells, n_spoints, n_sfacets = _cut_sub_triangulation_count_with_boundary(m,point_to_value)
  s, s_boundary = allocate_sub_triangulation_with_boundary(m,n_scells,n_spoints,n_sfacets)

  s, scell_to_inoutcut, s_boundary
end

function _cut_sub_triangulation_count(m,mpoint_to_value)
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

function _cut_sub_triangulation_count_with_boundary(m,mpoint_to_value)
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


