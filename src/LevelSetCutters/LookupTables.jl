## Table type and cache

"""
    struct LookupTable{D,T}

Marching-simplices table of a polytope `p` of dimension `D` with `N` vertices. Get it with
[`get_lookup_table`](@ref).

A case is the IN/OUT pattern of the level set at the vertices,
`case = 1 + Σ_{v OUT} 2^(v-1)` (see `compute_case`). Per case, the local points of
the cut polytope are its `N` vertices followed by one point per cut edge, in the order of
`case_to_cut_edges[case]`. The cutter appends the points of every cell in this same order, so
the local point ids of subcells and subfacets are valid in the cut cell.

Per-case fields (read only):
- `case_to_inoutcut`: `IN`, `OUT` or `CUT` (the whole polytope);
- `case_to_cut_edges`: the edges whose vertices have different sign, as pairs of local
  vertices, in the edge order of `get_faces(p,1,0)`;
- `case_to_subcell_to_points`: simplices covering `p` (positively oriented);
- `case_to_subcell_to_inout`: `IN` or `OUT` of every subcell;
- `case_to_subfacet_to_points`: simplices of the interface;
- `case_to_subfacet_to_orientation`: sign that makes the normal computed from the subfacet
  points point from IN to OUT.
"""
struct LookupTable{D,T}
  num_vertices::Int
  case_to_inoutcut::Vector{Int8}
  case_to_cut_edges::Vector{Vector{NTuple{2,Int}}}
  case_to_subcell_to_points::Vector{Vector{Vector{Int}}}
  case_to_subcell_to_inout::Vector{Vector{Int8}}
  case_to_subfacet_to_points::Vector{Vector{Vector{Int}}}
  case_to_subfacet_to_orientation::Vector{Vector{T}}
end

const LOOKUP_TABLES = Dict{Tuple{Vararg{Int}},LookupTable}()
const LOOKUP_TABLES_LOCK = ReentrantLock()

"""
    get_lookup_table(p::Polytope) -> LookupTable

The [`LookupTable`](@ref) of `p`, built on first use and cached for the session
(thread safe). Polytopes are identified by their extrusion (`Polytope` hashes by identity).
"""
function get_lookup_table(p::Polytope{D}) where D
  key = Tuple(get_extrusion(p).data)
  table = lock(LOOKUP_TABLES_LOCK) do
    get!(() -> LookupTable(p), LOOKUP_TABLES, key)
  end
  return table::LookupTable{D,Float64}
end

## Cases

isout(v) = v > 0

# The case of the values at the vertices of a cell: bit i-1 is set if vertex i is OUT
function compute_case(values)
  case = 1
  for (i,v) in enumerate(values)
    if isout(v)
      case += 2^(i-1)
    end
  end
  return case
end

function compute_case(cell_to_points::Table, point_to_value::AbstractVector, cell::Integer)
  compute_case(view(point_to_value,dataview(cell_to_points,cell)))
end

# IN if all values are ≤ 0, OUT if all are > 0, CUT otherwise (a zero value counts as IN)
function _find_in_out_or_cut(vertex_to_value)
  b = isout(first(vertex_to_value))
  for value in vertex_to_value
    if b != isout(value)
      return CUT
    end
  end
  return ifelse(b,  OUT,  IN)
end

# IN/OUT values at the vertices for a case (inverse of `compute_case`)
function _case_to_vertex_values(case,nvertices)
  [ isodd((case-1) >> (vertex-1)) ? OUT : IN for vertex in 1:nvertices ]
end

## Per-case queries

num_case_points(table::LookupTable,case) = table.num_vertices + length(table.case_to_cut_edges[case])

"""
    get_case_points(p::Polytope, case::Integer) -> Vector{<:Point}

Coordinates of the local points of `case` in the lookup table of `p`, as the table was
built: the vertices of `p`, then the midpoints of the cut edges. (The table does not store
them; the cutter places the cut points at the zero of the level set instead.)
"""
function get_case_points(p::Polytope, case::Integer)
  table = get_lookup_table(p)
  _case_points(get_vertex_coordinates(p),table.case_to_cut_edges[case])
end

## Building a table

function LookupTable(p::Polytope)
  vertex_to_coords = get_vertex_coordinates(p)
  edge_to_vertices = get_faces(p,1,0)
  D = num_dims(p)
  T = eltype(eltype(vertex_to_coords))
  nvertices = num_vertices(p)
  ncases = 2^nvertices

  case_to_inoutcut = Vector{Int8}(undef,ncases)
  case_to_cut_edges = Vector{Vector{NTuple{2,Int}}}(undef,ncases)
  case_to_subcell_to_points = Vector{Vector{Vector{Int}}}(undef,ncases)
  case_to_subcell_to_inout = Vector{Vector{Int8}}(undef,ncases)
  case_to_subfacet_to_points = Vector{Vector{Vector{Int}}}(undef,ncases)
  case_to_subfacet_to_orientation = Vector{Vector{T}}(undef,ncases)

  for case in 1:ncases
    vertex_to_value = _case_to_vertex_values(case,nvertices)
    cut_edges = _find_cut_edges(vertex_to_value,edge_to_vertices)
    point_to_coords = _case_points(vertex_to_coords,cut_edges)
    point_to_value = vcat(vertex_to_value,fill(INTERFACE,length(cut_edges)))
    subcell_to_points = _delaunay(point_to_coords)
    _ensure_positive_jacobians!(subcell_to_points,point_to_coords,p)
    subcell_to_inout = _compute_subcell_to_inout(subcell_to_points,point_to_value)
    subfacet_to_points, subfacet_to_normal = _find_subfacets(
      point_to_coords,subcell_to_points,subcell_to_inout,p)

    case_to_inoutcut[case] = _find_in_out_or_cut(vertex_to_value)
    case_to_cut_edges[case] = cut_edges
    case_to_subcell_to_points[case] = subcell_to_points
    case_to_subcell_to_inout[case] = subcell_to_inout
    case_to_subfacet_to_points[case] = subfacet_to_points
    case_to_subfacet_to_orientation[case] = _setup_subfacet_to_orientation(
      point_to_coords,subfacet_to_points,subfacet_to_normal)
  end

  LookupTable{D,T}(
    nvertices,
    case_to_inoutcut,
    case_to_cut_edges,
    case_to_subcell_to_points,
    case_to_subcell_to_inout,
    case_to_subfacet_to_points,
    case_to_subfacet_to_orientation)
end

function _find_cut_edges(vertex_to_value,edge_to_vertices)
  cut_edges = NTuple{2,Int}[]
  for (v1,v2) in edge_to_vertices
    if isout(vertex_to_value[v1]) != isout(vertex_to_value[v2])
      push!(cut_edges,(v1,v2))
    end
  end
  cut_edges
end

# The vertices, then the midpoints of the cut edges (the table is built on the midpoints;
# the cutter places the points at the interpolated zero of the level set)
function _case_points(vertex_to_coords,cut_edges)
  point_to_coords = copy(vertex_to_coords)
  for (v1,v2) in cut_edges
    push!(point_to_coords,0.5*(vertex_to_coords[v1]+vertex_to_coords[v2]))
  end
  point_to_coords
end

function _delaunay(points::Vector{Point{D,T}}) where {D,T}
  n = length(points)
  m = zeros(T,D,n)
  for (i,p) in enumerate(points)
    for (j,pj) in enumerate(p)
      m[j,i] = pj
    end
  end
  cells = delaunay(m)
  [ Vector{Int}(cells[:,k]) for k in 1:size(cells,2)]
end

function _ensure_positive_jacobians!(subcell_to_points,point_to_coords,p::Polytope)
  reffe = LagrangianRefFE(Float64,Simplex(p),1)
  shapefuns = get_shapefuns(reffe)
  vertex_to_coords = get_vertex_coordinates(p)
  shapefuns_grad = collect1d(evaluate(Broadcasting(∇)(shapefuns),vertex_to_coords)[1,:])
  _ensure_positive_jacobians_work!(subcell_to_points,point_to_coords,shapefuns_grad)
end

# Swaps the first two points of every cell with a negative jacobian
function _ensure_positive_jacobians_work!(cell_to_points,point_to_coords,shapefuns_grad)
  Ta = eltype(shapefuns_grad)
  Tb = eltype(point_to_coords)
  for cell in 1:length(cell_to_points)
    points = _cell_points(cell_to_points,cell)
    J = zero(outer(zero(Ta),zero(Tb)))
    for (i,point) in enumerate(points)
      J += outer(shapefuns_grad[i],point_to_coords[point])
    end
    if det(J) < 0
      points[1], points[2] = points[2], points[1]
    end
  end
end

# Mutable access to the points of a cell
_cell_points(cell_to_points::Table,cell) = dataview(cell_to_points,cell)
_cell_points(cell_to_points::AbstractVector,cell) = cell_to_points[cell]

# OUT if any point of the subcell is OUT (interface points are not)
function _compute_subcell_to_inout(subcell_to_points,point_to_value)
  Int8[ any(point -> isout(point_to_value[point]), points) ? OUT : IN
        for points in subcell_to_points ]
end

# The interface between IN and OUT subcells and its normals, via a Gridap model of the subcells
function _find_subfacets(point_to_coords,subcell_to_points,subcell_to_inout,p::Polytope)
  simplex = Simplex(p)
  reffe = LagrangianRefFE(Float64,simplex,1)
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
  q = get_cell_points(CellQuadrature(interface.⁺,0))
  n_q = collect(evaluate(get_normal_vector(interface.⁺),q))
  subfacet_to_normal = map(first,n_q)

  return subfacet_to_points, subfacet_to_normal
end

function _find_subfacet_to_points(
  subcell_to_points, subfacet_to_subcell, subfacet_to_lfacet, lfacet_to_lpoints)

  subfacet_to_points = Vector{Int}[]
  for subfacet in eachindex(subfacet_to_subcell)
    subcell = subfacet_to_subcell[subfacet]
    lfacet = subfacet_to_lfacet[subfacet]
    points = subcell_to_points[subcell]
    lpoints = lfacet_to_lpoints[lfacet]
    push!(subfacet_to_points,points[lpoints])
  end
  return subfacet_to_points
end

# Sign of the normal computed from the subfacet points w.r.t. the IN → OUT normal
function _setup_subfacet_to_orientation(point_to_coords,subfacet_to_points,subfacet_to_normal)
  T = eltype(eltype(point_to_coords))
  subfacet_to_orientation = zeros(T,length(subfacet_to_normal))
  for subfacet in eachindex(subfacet_to_points)
    v = _setup_normal(subfacet_to_points,point_to_coords,subfacet,0)
    subfacet_to_orientation[subfacet] = sign(v⋅subfacet_to_normal[subfacet])
  end
  return subfacet_to_orientation
end

## Normals of simplices (also used by the cutter and DifferentiableTriangulations)

# Unit normal of a subfacet from its points (shifted by `pointoffset`), not oriented
function _setup_normal(subfacet_to_points,point_to_coords::Vector{<:Point{1}},subfacet,pointoffset)
  T = eltype(eltype(point_to_coords))
  VectorValue(one(T))
end

function _setup_normal(
  subfacet_to_points,point_to_coords::Vector{<:Point{D}},subfacet,pointoffset) where D

  points = subfacet_to_points[subfacet]
  x1 = point_to_coords[points[1]+pointoffset]
  edges = ntuple(i -> point_to_coords[points[i+1]+pointoffset] - x1, Val(D-1))
  return _normal_vector(edges...)
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
