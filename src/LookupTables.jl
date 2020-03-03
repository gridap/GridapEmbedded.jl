
struct LookupTable{D,T}
  case_to_subcell_to_points::Vector{Vector{Vector{Int}}}
  case_to_subcell_to_inout::Vector{Vector{Int}}
  case_to_subfacet_to_subcell::Vector{Vector{Int}}
  case_to_subfacet_to_lfacet::Vector{Vector{Int}}
  case_to_subfacet_to_points::Vector{Vector{Vector{Int}}}
  case_to_subfacet_to_normal::Vector{Vector{VectorValue{D,T}}}
  case_to_point_to_coordinates::Vector{Vector{VectorValue{D,T}}} 
  case_to_inoutcut::Vector{Int}
end

function LookupTable(p::Polytope)
  _LookupTable(p)
end

#function is_cut_cell(table::LookupTable,ls::AbstractVector...)
#end

function num_cases(nvertices)
  2^nvertices
end

function isout(v)
  v > 0
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

const IN = -1
const OUT = 1
const INTERFACE = 0
const CUT = 0


# Helpers

function Simplex(p::Polytope)
  D = num_cell_dims(p)
  extrusion = tfill(TET_AXIS,Val{D}())
  ExtrusionPolytope(extrusion)
end

function Simplex(p::Polytope{2})
  TRI
end

function Simplex(p::Polytope{3})
  TET
end

function _compute_case(values)
  t = Table([collect(1:length(values)),])
  compute_case(t,values,1)
end

function _LookupTable(p::Polytope)

  nvertices = num_vertices(p)

  edge_to_vertices = get_faces(p,1,0)
  vertex_to_coords = get_vertex_coordinates(p)
  simplex = Simplex(p)
  reffe = LagrangianRefFE(Float64,simplex,1)
  D = num_cell_dims(simplex)
  lfacet_to_lpoints = get_faces(simplex,D-1,0)
  degree = 0

  V = eltype(get_vertex_coordinates(p))

  cis = CartesianIndices(Tuple(fill(2,nvertices)))

  ncases = length(cis)

  case_to_subcell_to_points = Vector{Vector{Vector{Int}}}(undef,ncases)
  case_to_subcell_to_inout = Vector{Vector{Int}}(undef,ncases)
  case_to_subfacet_to_subcell = Vector{Vector{Int}}(undef,ncases)
  case_to_subfacet_to_lfacet = Vector{Vector{Int}}(undef,ncases)
  case_to_subfacet_to_points = Vector{Vector{Vector{Int}}}(undef,ncases)
  case_to_subfacet_to_normal = Vector{Vector{V}}(undef,ncases)
  case_to_point_to_coordinates = Vector{Vector{V}}(undef,ncases)
  case_to_inoutcut = Vector{Int}(undef,ncases)

  for ci in cis

    vertex_to_value = _prepare_vertex_to_value(ci)
    point_to_coords, point_to_value = _compute_delaunay_points(vertex_to_value, vertex_to_coords, edge_to_vertices)
    subcell_to_points = delaunay(point_to_coords)
    subcell_to_inout = _compute_subcell_to_inout(subcell_to_points,point_to_value)
    subcell_to_ctype = fill(Int8(1),length(subcell_to_points))
    subgrid = UnstructuredGrid(point_to_coords,Table(subcell_to_points),[reffe],subcell_to_ctype)
    topo = GridTopology(subgrid)
    labels = FaceLabeling(topo)
    model = DiscreteModel(subgrid,topo,labels)
    interface = InterfaceTriangulation(model,collect(Bool,subcell_to_inout .== IN))
    subfacet_to_subcell = get_face_to_cell(interface.left)
    subfacet_to_lfacet = get_face_to_lface(interface.left)
    subfacet_to_points = _find_subfacet_to_points(
      subcell_to_points, subfacet_to_subcell, subfacet_to_lfacet, lfacet_to_lpoints)
    quad = CellQuadrature(interface,degree)
    q = get_coordinates(quad)
    n = get_normal_vector(interface)
    n_q = collect(evaluate(n,q))
    subfacet_to_normal = map(first,n_q)
    inoutcut = _find_in_out_or_cut(vertex_to_value)

    case = _compute_case(vertex_to_value)
    case_to_subcell_to_points[case] = subcell_to_points
    case_to_subcell_to_inout[case] = subcell_to_inout
    case_to_subfacet_to_subcell[case] = subfacet_to_subcell
    case_to_subfacet_to_lfacet[case] = subfacet_to_lfacet
    case_to_subfacet_to_points[case] = subfacet_to_points
    case_to_subfacet_to_normal[case] = subfacet_to_normal
    case_to_point_to_coordinates[case] = point_to_coords
    case_to_inoutcut[case] = inoutcut

  end

  LookupTable(
    case_to_subcell_to_points,
    case_to_subcell_to_inout,
    case_to_subfacet_to_subcell,
    case_to_subfacet_to_lfacet,
    case_to_subfacet_to_points,
    case_to_subfacet_to_normal,
    case_to_point_to_coordinates,
    case_to_inoutcut)
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
  subcell_to_points
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

function delaunay(points::Vector{Point{D,T}}) where {D,T}
  n = length(points)
  m = zeros(T,D,n)
  for (i,p) in enumerate(points)
    for (j,pj) in enumerate(p)
      m[j,i] = pj
    end
  end
  cells = delaunay(m)
  [ cells[:,k] for k in 1:size(cells,2)]
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

