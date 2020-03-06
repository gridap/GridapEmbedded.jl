
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

  a = cell_to_points.ptrs[cell]
  b = cell_to_points.ptrs[cell+1]-1
  all_in = true
  all_out = true
  for p in a:b
    point = cell_to_points.data[p]
    maxvalue = -1
    for point_to_value in ls_to_point_to_value
      value = point_to_value[point]
      maxvalue = max(maxvalue,value)
    end
    all_in = all_in && (!isout(maxvalue))
    all_out = all_out && isout(maxvalue)
  end
  if all_in
    return IN
  elseif all_out
    return OUT
  else
    return CUT
  end
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

struct SubTriangulation{Dc,Df,T}
  cell_table::LookupTable{Dc,T}
  facet_table::LookupTable{Df,T}
  cell_to_points::Table{Int,Int32}
  cell_to_inoutcut::Vector{Int8}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dc,T}}
  ls_to_point_to_value::Vector{Vector{T}}
end

struct FacetSubTriangulation{Dc,Df,T}
  facet_table::LookupTable{Df,T}
  facet_to_points::Table{Int,Int32}
  facet_to_normal::Vector{Point{Dc,T}}
  facet_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dc,T}}
  ls_to_point_to_value::Vector{Vector{T}}
end

function initial_sub_triangulation(grid::Grid,point_to_value::AbstractArray)
  initial_sub_triangulation(grid,[point_to_value,])
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
  cell_table = LookupTable(simplex)
  lpoint_to_lcoords = get_vertex_coordinates(p)
  _ensure_positive_jacobians!(ltcell_to_lpoints,lpoint_to_lcoords,cell_table.subcell_shapefuns_grad)
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
  facet_simplex = FacetSimplex(simplex)
  facet_table = LookupTable(facet_simplex)

  SubTriangulation(
    cell_table,
    facet_table,
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

function UnstructuredGrid(st::FacetSubTriangulation{Dc,Df}) where {Dc,Df}
  reffe = LagrangianRefFE(Float64,Simplex(Val{Df}()),1)
  cell_types = fill(Int8(1),length(st.facet_to_points))
  UnstructuredGrid(
    st.point_to_coords,
    st.facet_to_points,
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

function cut_sub_triangulation(st::SubTriangulation{Dc,Df,T}) where {Dc,Df,T}
  _st = st
  ls_to_fst = FacetSubTriangulation{Dc,Df,T}[]
  while length(_st.ls_to_point_to_value) > 0
    point_to_value = pop!(_st.ls_to_point_to_value)
    _st, _fst = _cut_sub_triangulation(_st,point_to_value)
    _fst = cut_sub_triangulation(_fst)
    push!(ls_to_fst,_fst)
  end
  _st, ls_to_fst
end

function cut_sub_triangulation(st::FacetSubTriangulation)
  _st = st
  while length(_st.ls_to_point_to_value) > 0
    point_to_value = pop!(_st.ls_to_point_to_value)
    _st = _cut_sub_triangulation(_st,point_to_value)
  end
  _st
end

function _cut_sub_triangulation(st::SubTriangulation,point_to_value)
  nrcells, nrpoints, nrfacets = _cut_sub_triangulation_count(st,point_to_value)
  rst, rfst = _allocate_new_sub_triangulation(st,nrcells,nrpoints,nrfacets)
  _fill_new_sub_triangulation!(rst,rfst,st,point_to_value)
  rst, rfst
end

function _cut_sub_triangulation(st::FacetSubTriangulation,point_to_value)
  nrfacets, nrpoints = _cut_sub_triangulation_count(st,point_to_value)
  rst = _allocate_new_sub_triangulation(st,nrfacets,nrpoints)
  _fill_new_sub_triangulation!(rst,st,point_to_value)
  rst
end

function _cut_sub_triangulation_count(st::SubTriangulation,point_to_value)
  nrcells = 0
  nrpoints = 0
  nrfacets = 0
  for cell in 1:length(st.cell_to_points)
    case = compute_case(st.cell_to_points,point_to_value,cell)
    nrcells += length(st.cell_table.case_to_subcell_to_inout[case])
    nrpoints += length(st.cell_table.case_to_point_to_coordinates[case])
    if st.cell_to_inoutcut[cell] != OUT
      nrfacets += length(st.cell_table.case_to_subfacet_to_subcell[case])
    end
  end
  nrcells, nrpoints, nrfacets
end

function _cut_sub_triangulation_count(st::FacetSubTriangulation,point_to_value)
  nrfacets = 0
  nrpoints = 0
  for facet in 1:length(st.facet_to_points)
    case = compute_case(st.facet_to_points,point_to_value,facet)
    if st.facet_table.case_to_inoutcut[case] != OUT
      nrfacets += st.facet_table.case_to_num_in[case]
      nrpoints += length(st.facet_table.case_to_point_to_coordinates[case])
    end
  end
  nrfacets, nrpoints
end

function _allocate_new_sub_triangulation(st::SubTriangulation{Dc,Df,T},nrcells,nrpoints,nrfacets) where {Dc,Df,T}

  nlp = st.cell_table.nlpoints
  rcell_to_rpoints_data = zeros(eltype(st.cell_to_points.data),nlp*nrcells)
  rcell_to_rpoints_ptrs = fill(eltype(st.cell_to_points.ptrs)(nlp),nrcells+1)
  length_to_ptrs!(rcell_to_rpoints_ptrs)
  rcell_to_rpoints = Table(rcell_to_rpoints_data,rcell_to_rpoints_ptrs)
  rcell_to_inoutcut = zeros(eltype(st.cell_to_inoutcut),nrcells)
  rcell_to_bgcell = zeros(eltype(st.cell_to_bgcell),nrcells)
  rpoint_to_coords = zeros(Point{Dc,T},nrpoints)
  ls_to_rpoint_to_value = [zeros(T,nrpoints) for i in 1:length(st.ls_to_point_to_value)]

  _st = SubTriangulation(
    st.cell_table,
    st.facet_table,
    rcell_to_rpoints,
    rcell_to_inoutcut,
    rcell_to_bgcell,
    rpoint_to_coords,
    ls_to_rpoint_to_value)

  nlpf = st.facet_table.nlpoints
  rfacet_to_rpoints_data = zeros(eltype(st.cell_to_points.data),nlpf*nrfacets)
  rfacet_to_rpoints_ptrs = fill(eltype(st.cell_to_points.ptrs)(nlpf),nrfacets+1)
  length_to_ptrs!(rfacet_to_rpoints_ptrs)
  rfacet_to_rpoints = Table(rfacet_to_rpoints_data,rfacet_to_rpoints_ptrs)
  rfacet_to_normal = zeros(VectorValue{Dc,T},nrfacets)
  rfacet_to_bgcell = zeros(eltype(st.cell_to_bgcell),nrfacets)

  _fst = FacetSubTriangulation(
    st.facet_table,
    rfacet_to_rpoints,
    rfacet_to_normal,
    rfacet_to_bgcell,
    rpoint_to_coords,
    [  rpoint_to_value for rpoint_to_value in ls_to_rpoint_to_value ])

  _st, _fst
end

function _allocate_new_sub_triangulation(st::FacetSubTriangulation{Dc,Df,T},nrfacets,nrpoints) where {Dc,Df,T}

  nlpf = st.facet_table.nlpoints
  rfacet_to_rpoints_data = zeros(eltype(st.facet_to_points.data),nlpf*nrfacets)
  rfacet_to_rpoints_ptrs = fill(eltype(st.facet_to_points.ptrs)(nlpf),nrfacets+1)
  length_to_ptrs!(rfacet_to_rpoints_ptrs)
  rfacet_to_rpoints = Table(rfacet_to_rpoints_data,rfacet_to_rpoints_ptrs)
  rfacet_to_normal = zeros(VectorValue{Dc,T},nrfacets)
  rfacet_to_bgcell = zeros(eltype(st.facet_to_bgcell),nrfacets)
  rpoint_to_coords = zeros(Point{Dc,T},nrpoints)
  ls_to_rpoint_to_value = [zeros(T,nrpoints) for i in 1:length(st.ls_to_point_to_value)]

  _fst = FacetSubTriangulation(
    st.facet_table,
    rfacet_to_rpoints,
    rfacet_to_normal,
    rfacet_to_bgcell,
    rpoint_to_coords,
    ls_to_rpoint_to_value)

  _fst
end

function _fill_new_sub_triangulation!(rst,rfst,st,point_to_value)
  rcell = 0
  rpoint = 0
  rfacet = 0
  q = 0
  z = 0
  for cell in 1:length(st.cell_to_points)

    pointoffset = rpoint
    a = st.cell_to_points.ptrs[cell]-1
    for lpoint in 1:st.cell_table.nlpoints
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

    if CUT == st.cell_table.case_to_inoutcut[case]
      for (ledge, lpoints) in enumerate(st.cell_table.ledge_to_lpoints)
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

    celloffset = rcell
    nsubcells = length(st.cell_table.case_to_subcell_to_inout[case])
    for subcell in 1:nsubcells
      rcell += 1

      if st.cell_to_inoutcut[cell] == OUT
        rst.cell_to_inoutcut[rcell] = OUT
      else
        rst.cell_to_inoutcut[rcell] = st.cell_table.case_to_subcell_to_inout[case][subcell]
      end

      rst.cell_to_bgcell[rcell] = st.cell_to_bgcell[cell]
      for subpoint in st.cell_table.case_to_subcell_to_points[case][subcell]
        q += 1
        rst.cell_to_points.data[q] = subpoint + pointoffset
      end
    end

    if st.cell_to_inoutcut[cell] != OUT
      nsubfacets = length(st.cell_table.case_to_subfacet_to_subcell[case])
      for subfacet in 1:nsubfacets
        rfacet += 1
        for subpoint in st.cell_table.case_to_subfacet_to_points[case][subfacet]
          z += 1
          rfst.facet_to_points.data[z] = subpoint + pointoffset
        end
        rfst.facet_to_bgcell[rfacet] = st.cell_to_bgcell[cell]
        normal = _setup_normal(
          st.cell_table.case_to_subfacet_to_points[case],
          rfst.point_to_coords,
          subfacet,pointoffset)
        orientation = st.cell_table.case_to_subfacet_to_orientation[case][subfacet]
        rfst.facet_to_normal[rfacet] = orientation*normal
      end
    end

  end
end

function _fill_new_sub_triangulation!(rfst::FacetSubTriangulation,st::FacetSubTriangulation,point_to_value)
  rpoint = 0
  rfacet = 0
  q = 0
  for facet in 1:length(st.facet_to_points)

    case = compute_case(st.facet_to_points,point_to_value,facet)
    if OUT == st.facet_table.case_to_inoutcut[case]
      continue
    end

    pointoffset = rpoint
    a = st.facet_to_points.ptrs[facet]-1
    for lpoint in 1:st.facet_table.nlpoints
      rpoint += 1
      point = st.facet_to_points.data[a+lpoint]
      coords = st.point_to_coords[point]
      rfst.point_to_coords[rpoint] = coords
      for (ils, point_to_val) in enumerate(st.ls_to_point_to_value)
        val = point_to_val[point]
        rfst.ls_to_point_to_value[ils][rpoint] = val
      end
    end


    if CUT == st.facet_table.case_to_inoutcut[case]
      for (ledge, lpoints) in enumerate(st.facet_table.ledge_to_lpoints)
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
          for (ils, point_to_val) in enumerate(st.ls_to_point_to_value)
            s1 = point_to_val[point1]
            s2 = point_to_val[point2]
            ds = s2-s1
            s = s1 + c1*ds
            rfst.ls_to_point_to_value[ils][rpoint] = s
          end
        end
      end
    end

    nsubfacets = length(st.facet_table.case_to_subcell_to_inout[case])
    for subfacet in 1:nsubfacets
      if st.facet_table.case_to_subcell_to_inout[case][subfacet] == IN
        rfacet += 1
        rfst.facet_to_bgcell[rfacet] = st.facet_to_bgcell[facet]
        rfst.facet_to_normal[rfacet] = st.facet_to_normal[facet]
        for subpoint in st.facet_table.case_to_subcell_to_points[case][subfacet]
          q += 1
          rfst.facet_to_points.data[q] = subpoint + pointoffset
        end
      end
    end

  end
end

function _setup_normal(subfacet_to_points,point_to_coords::Vector{<:Point{1}},subfacet,pointoffset)
  T = eltype(eltype(point_to_coords))
  VectorValue(one(T))
end

function _setup_normal(subfacet_to_points,point_to_coords::Vector{<:Point{2}},subfacet,pointoffset)
  subpoints = subfacet_to_points[subfacet]
  p1 = point_to_coords[subpoints[1]+pointoffset]
  p2 = point_to_coords[subpoints[2]+pointoffset]
  v1 = p2-p1
  _normal_vector(v1)
end

function _setup_normal(subfacet_to_points,point_to_coords::Vector{<:Point{3}},subfacet,pointoffset)
  subpoints = subfacet_to_points[subfacet]
  p1 = point_to_coords[subpoints[1]+pointoffset]
  p2 = point_to_coords[subpoints[2]+pointoffset]
  p3 = point_to_coords[subpoints[3]+pointoffset]
  v1 = p2-p1
  v2 = p3-p1
  _normal_vector(v1,v2)
end

function _setup_normal(subfacet_to_points,point_to_coords::Vector{<:Point{4}},subfacet,pointoffset)
  subpoints = subfacet_to_points[subfacet]
  p1 = point_to_coords[subpoints[1]+pointoffset]
  p2 = point_to_coords[subpoints[2]+pointoffset]
  p3 = point_to_coords[subpoints[3]+pointoffset]
  p4 = point_to_coords[subpoints[4]+pointoffset]
  v1 = p2-p1
  v2 = p3-p1
  v3 = p4-p1
  _normal_vector(v1,v2,v3)
end

function _normal_vector(u::VectorValue...)
  v = _orthogonal_vector(u...)
  m = sqrt(inner(v,v))
  if m < eps()
    return zero(v)
  else
    return v/m
  end
end

function _orthogonal_vector(v::VectorValue{2})
  w1 = v[2]
  w2 = -v[1]
  VectorValue(w1,w2)
end

function _orthogonal_vector(v1::VectorValue{3},v2::VectorValue{3})
  w1 = v1[2]*v2[3] - v1[3]*v2[2]
  w2 = v1[3]*v2[1] - v1[1]*v2[3]
  w3 = v1[1]*v2[2] - v1[2]*v2[1]
  VectorValue(w1,w2,w3)
end

function _orthogonal_vector(v1::VectorValue{4},v2::VectorValue{4},v3::VectorValue{4})
    v11 = v1[1]; v21 = v2[1]; v31 = v3[1]
    v12 = v1[2]; v22 = v2[2]; v32 = v3[2]
    v13 = v1[3]; v23 = v2[3]; v33 = v3[3]
    v14 = v1[4]; v24 = v2[4]; v34 = v3[4]
    w1 = (v12*v23*v34 - v12*v24*v33 - v13*v22*v34 + v13*v24*v32 + v14*v22*v33 - v14*v23*v32)
    w2 = (v11*v24*v33 - v11*v23*v34 + v13*v21*v34 - v13*v24*v31 - v14*v21*v33 + v14*v23*v31)
    w3 = (v11*v22*v34 - v11*v24*v32 - v12*v21*v34 + v12*v24*v31 + v14*v21*v32 - v14*v22*v31)
    w4 = (v11*v23*v32 - v11*v22*v33 + v12*v21*v33 - v12*v23*v31 - v13*v21*v32 + v13*v22*v31)
  VectorValue(w1,w2,w3,w4)
end

