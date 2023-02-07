
struct LookupTable{D,T}
  num_cases::Int
  case_to_subcell_to_points::Vector{Vector{Vector{Int}}}
  case_to_subcell_to_inout::Vector{Vector{Int}}
  case_to_num_subcells_in::Vector{Int}
  case_to_num_subcells_out::Vector{Int}
  case_to_subfacet_to_points::Vector{Vector{Vector{Int}}}
  case_to_subfacet_to_normal::Vector{Vector{VectorValue{D,T}}}
  case_to_subfacet_to_orientation::Vector{Vector{T}}
  case_to_point_to_coordinates::Vector{Vector{VectorValue{D,T}}}
  case_to_inoutcut::Vector{Int}
end

function LookupTable(p::Polytope)
  _build_lookup_table(p)
end

function compute_case(
  cell_to_points::Table, point_to_value::AbstractVector, cell::Integer)

  case = 1
  a = cell_to_points.ptrs[cell]
  b = cell_to_points.ptrs[cell+1]-1
  for (i,p) in enumerate(a:b)
    point = cell_to_points.data[p]
    v  = point_to_value[point]
    if  isout(v)
      case += 2^(i-1)
    end
  end
  case
end

function compute_case(values)
  t = Table([collect(1:length(values)),])
  compute_case(t,values,1)
end

function isout(v)
  v > 0
end

# Helpers

function _build_lookup_table(p::Polytope)

  table = _allocate_lookup_table(p)

  cis = _setup_cartesian_indices(p)

  for ci in cis

    vertex_to_value = _prepare_vertex_to_value(ci)
    case = compute_case(vertex_to_value)
    point_to_coords, point_to_value = _compute_delaunay_points(vertex_to_value,p)
    subcell_to_points = _delaunay(point_to_coords)
    _ensure_positive_jacobians!(subcell_to_points,point_to_coords,p)
    subcell_to_inout = _compute_subcell_to_inout(subcell_to_points,point_to_value)
    subfacet_to_points, subfacet_to_normal = _find_subfacets(
      point_to_coords,subcell_to_points,subcell_to_inout,p)
    subfacet_to_orientation = _setup_subfacet_to_orientation(
      point_to_coords,subfacet_to_points,subfacet_to_normal)
    inoutcut = _find_in_out_or_cut(vertex_to_value)

    table.case_to_subcell_to_points[case] = subcell_to_points
    table.case_to_subcell_to_inout[case] = subcell_to_inout
    table.case_to_num_subcells_in[case] = count(i-> i == IN, subcell_to_inout)
    table.case_to_num_subcells_out[case] = count(i-> i == OUT, subcell_to_inout)
    table.case_to_subfacet_to_points[case] = subfacet_to_points
    table.case_to_subfacet_to_normal[case] = subfacet_to_normal
    table.case_to_subfacet_to_orientation[case] = subfacet_to_orientation
    table.case_to_point_to_coordinates[case] = point_to_coords
    table.case_to_inoutcut[case] = inoutcut

  end

  table
end

function _allocate_lookup_table(p::Polytope)

  nvertices = num_vertices(p)
  ncases = 2^nvertices
  vertex_to_coords = get_vertex_coordinates(p)
  V = eltype(vertex_to_coords)
  D = length(V)
  T = eltype(V)

  case_to_subcell_to_points = Vector{Vector{Vector{Int}}}(undef,ncases)
  case_to_subcell_to_inout = Vector{Vector{Int}}(undef,ncases)
  case_to_num_subcells_in = Vector{Int}(undef,ncases)
  case_to_num_subcells_out = Vector{Int}(undef,ncases)
  case_to_subfacet_to_points = Vector{Vector{Vector{Int}}}(undef,ncases)
  case_to_subfacet_to_normal = Vector{Vector{VectorValue{D,T}}}(undef,ncases)
  case_to_subfacet_to_orientation = Vector{Vector{T}}(undef,ncases)
  case_to_point_to_coordinates = Vector{Vector{VectorValue{D,T}}}(undef,ncases)
  case_to_inoutcut = Vector{Int}(undef,ncases)

  LookupTable(
    ncases,
    case_to_subcell_to_points,
    case_to_subcell_to_inout,
    case_to_num_subcells_in,
    case_to_num_subcells_out,
    case_to_subfacet_to_points,
    case_to_subfacet_to_normal,
    case_to_subfacet_to_orientation,
    case_to_point_to_coordinates,
    case_to_inoutcut)

end

function _setup_cartesian_indices(p::Polytope)
  nvertices = num_vertices(p)
  cis = CartesianIndices(Tuple(fill(2,nvertices)))
  cis
end

function _prepare_vertex_to_value(ci)
  nvertices = length(ci)
  vertex_to_value = fill(IN,nvertices)
  for (vertex,c) in enumerate(Tuple(ci))
    if c == 1
       vertex_to_value[vertex] = IN
    else
       vertex_to_value[vertex] = OUT
    end
  end
  vertex_to_value
end

function _compute_delaunay_points(vertex_to_value,p::Polytope)
  edge_to_vertices = get_faces(p,1,0)
  vertex_to_coords = get_vertex_coordinates(p)
  _compute_delaunay_points(vertex_to_value, vertex_to_coords, edge_to_vertices)
end

function _compute_delaunay_points(vertex_to_value, vertex_to_coords, edge_to_vertices)
  point_to_coords = copy(vertex_to_coords)
  point_to_value = copy(vertex_to_value)
  for vertices in edge_to_vertices
    v1 = vertex_to_value[vertices[1]]
    v2 = vertex_to_value[vertices[2]]
    if isout(v1) != isout(v2)
      p1 = vertex_to_coords[vertices[1]]
      p2 = vertex_to_coords[vertices[2]]
      p = 0.5*(p1+p2)
      push!(point_to_coords,p)
      push!(point_to_value,INTERFACE)
    end
  end
  point_to_coords, point_to_value
end

function _delaunay(points::Vector{Point{D,T}}) where {D,T}
  n = length(points)
  m = zeros(T,n,D)
  for (i,p) in enumerate(points)
    for (j,pj) in enumerate(p)
      m[i,j] = pj
    end
  end
  cells = delaunay(m).simplices
  [ Vector{Int}(cells[k,:]) for k in 1:size(cells,1)]
end

function _ensure_positive_jacobians!(subcell_to_points,point_to_coords,p::Polytope)

  simplex = Simplex(p)
  order = 1
  reffe = LagrangianRefFE(Float64,simplex,order)
  D = num_cell_dims(simplex)
  lfacet_to_lpoints = get_faces(simplex,D-1,0)
  shapefuns = get_shapefuns(reffe)
  vertex_to_coords = get_vertex_coordinates(p)
  shapefuns_grad = collect1d(evaluate(Broadcasting(∇)(shapefuns),vertex_to_coords)[1,:])
   _ensure_positive_jacobians_work!(subcell_to_points,point_to_coords,shapefuns_grad)
end

function _ensure_positive_jacobians_work!(subcell_to_points,point_to_coords,shapefuns_grad)
  for (subcell,points) in enumerate(subcell_to_points)
    Ta = eltype(shapefuns_grad)
    Tb = eltype(point_to_coords)
    J = zero(outer(zero(Ta),zero(Tb)))
    for (i,point) in enumerate(points)
      J += outer(shapefuns_grad[i],point_to_coords[point])
    end
    dV = det(J)
    if dV < 0
      n1 = subcell_to_points[subcell][1]
      n2 = subcell_to_points[subcell][2]
      subcell_to_points[subcell][1] = n2
      subcell_to_points[subcell][2] = n1
    end
  end
end

function _ensure_positive_jacobians_work!(subcell_to_points::Table,point_to_coords,shapefuns_grad)
  Ta = eltype(shapefuns_grad)
  Tb = eltype(point_to_coords)
  for (subcell,points) in enumerate(subcell_to_points)
    J = zero(outer(zero(Ta),zero(Tb)))
    for (i,point) in enumerate(points)
      J += outer(shapefuns_grad[i],point_to_coords[point])
    end
    dV = det(J)
    if dV < 0
      a = subcell_to_points.ptrs[subcell]
      n1 = subcell_to_points.data[a]
      n2 = subcell_to_points.data[a+1]
      subcell_to_points.data[a] = n2
      subcell_to_points.data[a+1] = n1
    end
  end
end

function  _compute_subcell_to_inout(subcell_to_points,point_to_value)
  nsubcells = length(subcell_to_points)
  subcell_to_inout = fill(IN,nsubcells)
  for subcell in 1:nsubcells
    points = subcell_to_points[subcell]
    for point in points
      value = point_to_value[point]
      if isout(value)
        subcell_to_inout[subcell] = OUT
        break
      end
    end
  end
  subcell_to_inout
end

function _find_subfacets(point_to_coords,subcell_to_points, subcell_to_inout,p::Polytope)

  simplex = Simplex(p)
  order = 1
  reffe = LagrangianRefFE(Float64,simplex,order)
  D = num_cell_dims(simplex)
  lfacet_to_lpoints = get_faces(simplex,D-1,0)

  subcell_to_ctype = fill(Int8(1),length(subcell_to_points))
  subgrid = UnstructuredGrid(point_to_coords,Table(subcell_to_points),[reffe],subcell_to_ctype)
  topo = GridTopology(subgrid)
  labels = FaceLabeling(topo)
  model = DiscreteModel(subgrid,topo,labels)
  interface = InterfaceTriangulation(model,collect(Bool,subcell_to_inout .== IN))
  subfacet_to_subcell = interface.⁺.glue.face_to_cell
  subfacet_to_lfacet = interface.⁺.glue.face_to_lface
  subfacet_to_points = _find_subfacet_to_points(
    subcell_to_points, subfacet_to_subcell, subfacet_to_lfacet, lfacet_to_lpoints)
  degree = 0
  quad = CellQuadrature(interface.⁺,degree)
  q = get_cell_points(quad)
  n = get_normal_vector(interface.⁺)
  n_q = collect(evaluate(n,q))
  subfacet_to_normal = map(first,n_q)

  subfacet_to_points, subfacet_to_normal
end

function _find_subfacet_to_points(
  subcell_to_points,subfacet_to_subcell,subfacet_to_lfacet, lfacet_to_lpoints)

  subfacet_to_points = Vector{Int}[]
  for subfacet in 1:length(subfacet_to_subcell)
    subcell = subfacet_to_subcell[subfacet]
    lfacet = subfacet_to_lfacet[subfacet]
    points = subcell_to_points[subcell]
    lpoints = lfacet_to_lpoints[lfacet]
    push!(subfacet_to_points,points[lpoints])
  end
  subfacet_to_points
end

function _setup_subfacet_to_orientation(point_to_coords,subfacet_to_points,subfacet_to_normal)
  T = eltype(eltype(point_to_coords))
  subfacet_to_orientation = zeros(T,length(subfacet_to_normal))
  for (subfacet, points) in enumerate(subfacet_to_points)
    v = _setup_normal(subfacet_to_points,point_to_coords,subfacet,0)
    n = subfacet_to_normal[subfacet]
    a = v⋅n
    orientation = sign(a)
    subfacet_to_orientation[subfacet] = orientation
  end
  subfacet_to_orientation
end

function _find_in_out_or_cut(vertex_to_value)
  b = isout(first(vertex_to_value))
  for value in vertex_to_value
    if b != isout(value)
      return CUT
    end
  end
  if b
    return OUT
  else
    return IN
  end
  return -1
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
