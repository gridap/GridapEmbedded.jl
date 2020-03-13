
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

function UnstructuredGrid(st::FacetSubTriangulation{Dc}) where Dc
  reffe = LagrangianRefFE(Float64,Simplex(Val{Dc-1}()),1)
  cell_types = fill(Int8(1),length(st.facet_to_points))
  UnstructuredGrid(
    st.point_to_coords,
    st.facet_to_points,
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

function writevtk(st::FacetSubTriangulation,filename::String)
  ug = UnstructuredGrid(st)
  degree = 0
  quad = CellQuadrature(ug,degree)
  dS = integrate(1,ug,quad)
  write_vtk_file(ug,filename,celldata=[
    "normal"=>st.facet_to_normal,
    "bgcell"=>st.facet_to_bgcell,
    "dS"=>dS])
end

function initial_sub_triangulation(grid::Grid,geom::AnalyticalGeometry)
  geom_x = discretize(geom,grid)
  initial_sub_triangulation(grid,geom_x)
end

function initial_sub_triangulation(grid::Grid,geom::DiscreteGeometry)
  ugrid = UnstructuredGrid(grid)
  _initial_sub_triangulation(ugrid,geom)
end

# Helpers

struct MiniCell
  num_points::Int
  edge_to_points::Vector{Vector{Int}}
end

function cut_sub_triangulation(st::SubTriangulation{Dc,T},geom) where {Dc,T}
  _st = st
  ls_to_point_to_value_copy = [ point_to_value for point_to_value in geom.ls_to_point_to_value ]
  _geom = DiscreteGeometry(geom,ls_to_point_to_value_copy)

  ls_to_fst = FacetSubTriangulation{Dc,T}[]
  while length(_geom.ls_to_point_to_value) > 0
    point_to_value = pop!(_geom.ls_to_point_to_value)
    _st, _geom, _fst, _fgeom = _cut_sub_triangulation(_st,_geom,point_to_value)
    _fst = cut_sub_triangulation(_fst,_fgeom)
    push!(ls_to_fst,_fst)
  end
  _st, ls_to_fst
end

function cut_sub_triangulation(st::FacetSubTriangulation,geom)
  _st = st
  ls_to_point_to_value_copy = [ point_to_value for point_to_value in geom.ls_to_point_to_value ]
  _geom = DiscreteGeometry(geom,ls_to_point_to_value_copy)
  while length(_geom.ls_to_point_to_value) > 0
    point_to_value = pop!(_geom.ls_to_point_to_value)
    _st, _geom = _cut_sub_triangulation(_st,_geom,point_to_value)
  end
  _st
end

function _cut_sub_triangulation(st::SubTriangulation,geom,point_to_value)
  refcell, reffacet, cell_table, facet_table = _setup_tables(st)
  nrcells, nrpoints, nrfacets = _cut_sub_triangulation_count(st,geom,cell_table,point_to_value)
  rst, rstgeom, rfst, rfstgeom = _allocate_new_sub_triangulation(st,geom,refcell,reffacet,nrcells,nrpoints,nrfacets)
  _fill_new_sub_triangulation!(rst,rstgeom,rfst,cell_table,refcell,st,geom,point_to_value)
  rst, rstgeom, rfst, rfstgeom
end

function _cut_sub_triangulation(st::FacetSubTriangulation,geom,point_to_value)
  reffacet, facet_table = _setup_tables(st)
  nrfacets, nrpoints = _cut_sub_triangulation_count(st,geom,facet_table,point_to_value)
  rst, rstgeom = _allocate_new_sub_triangulation(st,geom,reffacet,nrfacets,nrpoints)
  _fill_new_sub_triangulation!(rst,rstgeom,facet_table,reffacet,st,geom,point_to_value)
  rst, rstgeom
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

function _setup_tables(st::FacetSubTriangulation{D}) where D
  pf = Simplex(Val{D-1}())
  facet = MiniCell(num_vertices(pf),get_faces(pf,1,0))
  facet_table = LookupTable(pf)
  facet, facet_table
end

function _cut_sub_triangulation_count(st::SubTriangulation,geom,cell_table,point_to_value)
  nrcells = 0
  nrpoints = 0
  nrfacets = 0
  if geom.intersection
    LOC = OUT
  else
    LOC = IN
  end
  for cell in 1:length(st.cell_to_points)
    case = compute_case(st.cell_to_points,point_to_value,cell)
    nrcells += length(cell_table.case_to_subcell_to_inout[case])
    nrpoints += length(cell_table.case_to_point_to_coordinates[case])
    if st.cell_to_inoutcut[cell] != LOC
      nrfacets += length(cell_table.case_to_subfacet_to_normal[case])
    end
  end
  nrcells, nrpoints, nrfacets
end

function _cut_sub_triangulation_count(st::FacetSubTriangulation,geom,facet_table,point_to_value)
  nrfacets = 0
  nrpoints = 0
  for facet in 1:length(st.facet_to_points)
    case = compute_case(st.facet_to_points,point_to_value,facet)
    if geom.intersection
      if facet_table.case_to_inoutcut[case] != OUT
        nrfacets += facet_table.case_to_num_subcells_in[case]
        nrpoints += length(facet_table.case_to_point_to_coordinates[case])
      end
    else
      if facet_table.case_to_inoutcut[case] != IN
        nrfacets += facet_table.case_to_num_subcells_out[case]
        nrpoints += length(facet_table.case_to_point_to_coordinates[case])
      end
    end
  end
  nrfacets, nrpoints
end

function _allocate_new_sub_triangulation(
  st::SubTriangulation{Dc,T},geom,refcell,reffacet,nrcells,nrpoints,nrfacets) where {Dc,T}

  nlp = refcell.num_points
  rcell_to_rpoints_data = zeros(eltype(st.cell_to_points.data),nlp*nrcells)
  rcell_to_rpoints_ptrs = fill(eltype(st.cell_to_points.ptrs)(nlp),nrcells+1)
  length_to_ptrs!(rcell_to_rpoints_ptrs)
  rcell_to_rpoints = Table(rcell_to_rpoints_data,rcell_to_rpoints_ptrs)
  rcell_to_inoutcut = zeros(eltype(st.cell_to_inoutcut),nrcells)
  rcell_to_bgcell = zeros(eltype(st.cell_to_bgcell),nrcells)
  rpoint_to_coords = zeros(Point{Dc,T},nrpoints)
  rpoint_to_rcoords = zeros(Point{Dc,T},nrpoints)
  ls_to_rpoint_to_value = [zeros(T,nrpoints) for i in 1:length(geom.ls_to_point_to_value)]

  _st = SubTriangulation(
    rcell_to_rpoints,
    rcell_to_inoutcut,
    rcell_to_bgcell,
    rpoint_to_coords,
    rpoint_to_rcoords)

  _st_geom = DiscreteGeometry(geom,ls_to_rpoint_to_value)

  nlpf = reffacet.num_points
  rfacet_to_rpoints_data = zeros(eltype(st.cell_to_points.data),nlpf*nrfacets)
  rfacet_to_rpoints_ptrs = fill(eltype(st.cell_to_points.ptrs)(nlpf),nrfacets+1)
  length_to_ptrs!(rfacet_to_rpoints_ptrs)
  rfacet_to_rpoints = Table(rfacet_to_rpoints_data,rfacet_to_rpoints_ptrs)
  rfacet_to_normal = zeros(VectorValue{Dc,T},nrfacets)
  rfacet_to_bgcell = zeros(eltype(st.cell_to_bgcell),nrfacets)

  ls_to_rpoint_to_value_copy = [ rpoint_to_value for rpoint_to_value in ls_to_rpoint_to_value ]

  _fst = FacetSubTriangulation(
    rfacet_to_rpoints,
    rfacet_to_normal,
    rfacet_to_bgcell,
    rpoint_to_coords,
    rpoint_to_rcoords)

  _fst_geom = DiscreteGeometry(geom,ls_to_rpoint_to_value_copy)

  _st, _st_geom, _fst, _fst_geom
end

function _allocate_new_sub_triangulation(
  st::FacetSubTriangulation{Dc,T},geom,reffacet,nrfacets,nrpoints) where {Dc,T}

  nlpf = reffacet.num_points
  rfacet_to_rpoints_data = zeros(eltype(st.facet_to_points.data),nlpf*nrfacets)
  rfacet_to_rpoints_ptrs = fill(eltype(st.facet_to_points.ptrs)(nlpf),nrfacets+1)
  length_to_ptrs!(rfacet_to_rpoints_ptrs)
  rfacet_to_rpoints = Table(rfacet_to_rpoints_data,rfacet_to_rpoints_ptrs)
  rfacet_to_normal = zeros(VectorValue{Dc,T},nrfacets)
  rfacet_to_bgcell = zeros(eltype(st.facet_to_bgcell),nrfacets)
  rpoint_to_coords = zeros(Point{Dc,T},nrpoints)
  rpoint_to_rcoords = zeros(Point{Dc,T},nrpoints)
  ls_to_rpoint_to_value = [zeros(T,nrpoints) for i in 1:length(geom.ls_to_point_to_value)]

  _fst = FacetSubTriangulation(
    rfacet_to_rpoints,
    rfacet_to_normal,
    rfacet_to_bgcell,
    rpoint_to_coords,
    rpoint_to_rcoords)

  _fst_geom = DiscreteGeometry(geom,ls_to_rpoint_to_value)

  _fst, _fst_geom
end

function _fill_new_sub_triangulation!(
  rst::SubTriangulation,rstgeom,rfst,cell_table,refcell,st,geom,point_to_value)

  rcell = 0
  rpoint = 0
  rfacet = 0
  q = 0
  z = 0
  if geom.intersection
    LOC = OUT
  else
    LOC = IN
  end
  for cell in 1:length(st.cell_to_points)

    pointoffset = rpoint
    a = st.cell_to_points.ptrs[cell]-1
    for lpoint in 1:refcell.num_points
      rpoint += 1
      point = st.cell_to_points.data[a+lpoint]
      coords = st.point_to_coords[point]
      rst.point_to_coords[rpoint] = coords
      for (ils, point_to_val) in enumerate(geom.ls_to_point_to_value)
        val = point_to_val[point]
        rstgeom.ls_to_point_to_value[ils][rpoint] = val
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
          for (ils, point_to_val) in enumerate(geom.ls_to_point_to_value)
            s1 = point_to_val[point1]
            s2 = point_to_val[point2]
            ds = s2-s1
            s = s1 + c1*ds
            rstgeom.ls_to_point_to_value[ils][rpoint] = s
          end
        end
      end
    end

    celloffset = rcell
    nsubcells = length(cell_table.case_to_subcell_to_inout[case])
    for subcell in 1:nsubcells
      rcell += 1

      if st.cell_to_inoutcut[cell] == LOC
        rst.cell_to_inoutcut[rcell] = LOC
      else
        rst.cell_to_inoutcut[rcell] = cell_table.case_to_subcell_to_inout[case][subcell]
      end

      rst.cell_to_bgcell[rcell] = st.cell_to_bgcell[cell]
      for subpoint in cell_table.case_to_subcell_to_points[case][subcell]
        q += 1
        rst.cell_to_points.data[q] = subpoint + pointoffset
      end
    end

    if st.cell_to_inoutcut[cell] != LOC
      nsubfacets = length(cell_table.case_to_subfacet_to_normal[case])
      for subfacet in 1:nsubfacets
        rfacet += 1
        for subpoint in cell_table.case_to_subfacet_to_points[case][subfacet]
          z += 1
          rfst.facet_to_points.data[z] = subpoint + pointoffset
        end
        rfst.facet_to_bgcell[rfacet] = st.cell_to_bgcell[cell]
        normal = _setup_normal(
          cell_table.case_to_subfacet_to_points[case],
          rfst.point_to_coords,
          subfacet,pointoffset)
        orientation = cell_table.case_to_subfacet_to_orientation[case][subfacet]
        rfst.facet_to_normal[rfacet] = orientation*normal
      end
    end

  end
end

function _fill_new_sub_triangulation!(
  rfst::FacetSubTriangulation,rfstgeom,facet_table,reffacet,st::FacetSubTriangulation,geom,point_to_value)
  rpoint = 0
  rfacet = 0
  q = 0
  if geom.intersection
    LOC = OUT
  else
    LOC = IN
  end
  for facet in 1:length(st.facet_to_points)

    case = compute_case(st.facet_to_points,point_to_value,facet)
    if LOC == facet_table.case_to_inoutcut[case]
      continue
    end

    pointoffset = rpoint
    a = st.facet_to_points.ptrs[facet]-1
    for lpoint in 1:reffacet.num_points
      rpoint += 1
      point = st.facet_to_points.data[a+lpoint]
      coords = st.point_to_coords[point]
      rfst.point_to_coords[rpoint] = coords
      for (ils, point_to_val) in enumerate(geom.ls_to_point_to_value)
        val = point_to_val[point]
        rfstgeom.ls_to_point_to_value[ils][rpoint] = val
      end
    end


    if CUT == facet_table.case_to_inoutcut[case]
      for (ledge, lpoints) in enumerate(reffacet.edge_to_points)
        point1 = st.facet_to_points.data[a+lpoints[1]]
        point2 = st.facet_to_points.data[a+lpoints[2]]
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
          rfst.point_to_coords[rpoint] = p
          for (ils, point_to_val) in enumerate(geom.ls_to_point_to_value)
            s1 = point_to_val[point1]
            s2 = point_to_val[point2]
            ds = s2-s1
            s = s1 + c1*ds
            rfstgeom.ls_to_point_to_value[ils][rpoint] = s
          end
        end
      end
    end

    nsubfacets = length(facet_table.case_to_subcell_to_inout[case])
    for subfacet in 1:nsubfacets
      if facet_table.case_to_subcell_to_inout[case][subfacet] != LOC
        rfacet += 1
        rfst.facet_to_bgcell[rfacet] = st.facet_to_bgcell[facet]
        rfst.facet_to_normal[rfacet] = st.facet_to_normal[facet]
        for subpoint in facet_table.case_to_subcell_to_points[case][subfacet]
          q += 1
          rfst.facet_to_points.data[q] = subpoint + pointoffset
        end
      end
    end

  end
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
  cutgeom = DiscreteGeometry(geom,ls_to_point_to_value)

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

  subgeom = DiscreteGeometry(cutgeom, ls_to_tpoint_to_value)

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
