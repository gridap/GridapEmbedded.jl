# Copied from background-gmg, src/cutters/cellcutter/CutFaces.jl (the recursive cell cutter, originally
# PolytopalCutters.jl). Changes are marked with [changed].

## The faces of a cut cell
#
# The leaves of a `CutCell` are the cells of the arrangement of the level-set planes inside
# the cell, so together with their faces they form a face-to-face complex. A d-face of this
# complex lies on exactly D-d planes (non-degeneracy) and is identified by them and by its
# phase: `IN`/`OUT` for the level sets it is not on, `INTERFACE` for those it is on (the
# region of the flat of its planes with a given sign vector is convex, so it is one face).
#
# The d-faces are read off the vertices of the leaves: in a simple polytope every set of D-d of
# the D planes of a vertex is a d-face through it, made of the vertices on all of them. Every
# vertex emits these keys (planes, phase) and a dictionary groups the vertices by key, so a
# face shared by two leaves (on a level set) comes out once. Faces can be selected by
# dimension and by the planes they lie on (`CellFaceExtractor`).
#
# Every face is split into d-simplices by a pulling triangulation: cone from its smallest
# vertex id over the triangulations of its facets not containing it. A face is triangulated
# from its vertex set only, so neighbouring faces of the cell are triangulated conformingly.

"""
    struct CutFaces{N}

The d-faces of a cut cell for one dimension d (an entry of `CutCellFaces`):

- `face_planes`: the `D-d` sorted plane ids of every face, flat (stride `D-d`);
- `face_phase`: `IN`/`OUT` per level set, `INTERFACE` for the level sets the face is on;
- `face_ptrs`, `face_vertices`: the vertices of every face (a table);
- `simplex_vertices`: the `d+1` vertices of every simplex, flat (stride `d+1`);
- `face_simplex_ptrs`: the simplices of face `f` are `face_simplex_ptrs[f]:face_simplex_ptrs[f+1]-1`.

Vertex ids are those of the `CutCell`.
"""
struct CutFaces{N}
  d::Int
  nplanes::Int
  face_planes::Vector{Int32}
  face_phase::Vector{NTuple{N,Int8}}
  face_ptrs::Vector{Int32}
  face_vertices::Vector{Int32}
  simplex_vertices::Vector{Int32}
  face_simplex_ptrs::Vector{Int32}
end

CutFaces{N}(d::Integer,D::Integer) where N =
  CutFaces{N}(d,D-d,Int32[],NTuple{N,Int8}[],Int32[1],Int32[],Int32[],Int32[1])

function _reset!(f::CutFaces)
  empty!(f.face_planes)
  empty!(f.face_phase)
  resize!(f.face_ptrs,1)
  empty!(f.face_vertices)
  empty!(f.simplex_vertices)
  resize!(f.face_simplex_ptrs,1)
  return f
end

"""
    struct CutCellFaces{D,N,T}

The faces of every dimension of a `CutCell` (`get_cut_cell`, the last one evaluated), with
their simplices (vertex coordinates: those of the cut cell). It is the cache and the output of
`CellFaceExtractor`, so it is overwritten for every cell. Accessors: `num_faces`,
`get_face_planes`, `get_face_phase`, `get_face_vertices`, `get_face_normal`,
`num_simplices`, `get_simplex_vertices`, `get_face_simplices`, `get_face_quadrature`.
"""
mutable struct CutCellFaces{D,N,T,C<:CutCell{D,N,T}}
  cell::C
  faces::Vector{CutFaces{N}}
  # Work buffers of the extraction: face of every key (planes padded with 0, phase), and the
  # (face, vertex) incidences
  face_index::Dict{Tuple{NTuple{D,Int32},NTuple{N,Int8}},Int32}
  face_vertex_pairs::Vector{Tuple{Int32,Int32}}
  # Work buffers of the triangulation: vertices per depth, planes and apices (stacks)
  vertex_buffers::Vector{Vector{Int32}}
  plane_stack::Vector{Int32}
  apex_stack::Vector{Int32}
end

function CutCellFaces(c::CutCell{D,N,T}) where {D,N,T}
  CutCellFaces(c,[CutFaces{N}(d,D) for d in 0:D],
    Dict{Tuple{NTuple{D,Int32},NTuple{N,Int8}},Int32}(),Tuple{Int32,Int32}[],
    [Int32[] for _ in 0:D],Int32[],Int32[])
end

"""
    CellFaceExtractor(; dims=nothing, planes=nothing)

`Map` from a `CutCell` to its `CutCellFaces`:

    cell_faces = lazy_map(CellFaceExtractor(),cell_cuts)

Only the faces of dimension in `dims` lying on planes in `planes` (plane ids, e.g. `1:F` for
the faces on the facets of the cell) are extracted; `nothing`: all. The leaves (d = D) lie on
no plane.
"""
struct CellFaceExtractor{A,B} <: Map
  dims::A
  planes::B
end
CellFaceExtractor(; dims=nothing, planes=nothing) = CellFaceExtractor(dims,planes)

_selected(::Nothing,x) = true
_selected(xs,x) = x in xs

Arrays.return_cache(::CellFaceExtractor,c::CutCell) = CutCellFaces(c)

function Arrays.evaluate!(cf::CutCellFaces{D,N,T,C},k::CellFaceExtractor,c::C) where {D,N,T,C}
  cf.cell = c
  _extract!(cf,k)
end

## Accessors

get_cut_cell(cf::CutCellFaces) = cf.cell

get_vertex_coordinates(cf::CutCellFaces,v::Integer) = get_vertex_coordinates(cf.cell,v)

num_faces(cf::CutCellFaces,d::Integer) = length(cf.faces[d+1].face_phase)

function get_face_planes(cf::CutCellFaces,d::Integer,f::Integer)
  fs = cf.faces[d+1]
  k = fs.nplanes
  view(fs.face_planes,(f-1)*k+1:f*k)
end

get_face_phase(cf::CutCellFaces,d::Integer,f::Integer) = cf.faces[d+1].face_phase[f]

function get_face_vertices(cf::CutCellFaces,d::Integer,f::Integer)
  fs = cf.faces[d+1]
  view(fs.face_vertices,fs.face_ptrs[f]:fs.face_ptrs[f+1]-1)
end

num_simplices(cf::CutCellFaces,d::Integer) = length(cf.faces[d+1].simplex_vertices) ÷ (d+1)

get_simplex_vertices(cf::CutCellFaces,d::Integer,s::Integer) =
  view(cf.faces[d+1].simplex_vertices,(s-1)*(d+1)+1:s*(d+1))

function get_face_simplices(cf::CutCellFaces,d::Integer,f::Integer)
  ptrs = cf.faces[d+1].face_simplex_ptrs
  ptrs[f]:ptrs[f+1]-1
end

"""
    get_face_normal(cf::CutCellFaces{D},f::Integer)

Unit normal (reference coordinates) of the (D-1)-face `f`: the normal of its plane, towards
`φ > 0` on a level set and outwards on a cell facet.
"""
get_face_normal(cf::CutCellFaces{D},f::Integer) where D = get_plane_normal(cf.cell,only(get_face_planes(cf,D-1,f)))

"""
    get_face_quadrature(cf::CutCellFaces,d,f,degree) -> (points, weights)

Gauss rule of `degree` on the d-face `f`, in reference coordinates of the cell: Gridap's
rule on the reference simplex mapped onto every simplex of the face (weights: d-measure).
"""
function get_face_quadrature(cf::CutCellFaces{D,N,T},d::Integer,f::Integer,degree::Integer) where {D,N,T}
  P = VectorValue{D,T}
  xs, ws = P[], T[]
  if d == 0
    push!(xs,get_vertex_coordinates(cf,only(get_face_vertices(cf,0,f))))
    push!(ws,one(T))
    return xs, ws
  end
  q = Quadrature((SEGMENT,TRI,TET)[d],degree)
  for s in get_face_simplices(cf,d,f)
    X = map(v -> get_vertex_coordinates(cf,v),get_simplex_vertices(cf,d,s))
    J = [(X[j+1]-X[1])[i] for i in 1:D, j in 1:d]
    μ = abs(det(qr(J).R)) # = √det(JᵀJ), accurate also for (shifted) zero-measure simplices
    for (ξ,w) in zip(get_coordinates(q),get_weights(q))
      push!(xs,X[1] + sum(j -> ξ[j]*(X[j+1]-X[1]),1:d))
      push!(ws,w*μ)
    end
  end
  xs, ws
end

## Extraction

function _extract!(cf::CutCellFaces{D},k::CellFaceExtractor) where D
  c = cf.cell
  for d in 0:D
    fs = _reset!(cf.faces[d+1])
    _selected(k.dims,d) || continue
    if d == D
      foreach(n -> _add_leaf!(fs,c,n),get_leaves(c))
    else
      _add_faces!(cf,fs,k.planes)
    end
    for f in 1:length(fs.face_phase)
      _triangulate_face!(cf,fs,f)
    end
  end
  return cf
end

function _add_leaf!(fs::CutFaces,c::CutCell,n::Integer)
  append!(fs.face_vertices,get_node_vertices(c,n))
  push!(fs.face_ptrs,length(fs.face_vertices)+1)
  push!(fs.face_phase,get_node_phase(c,n))
end

# The d-faces of the leaves, d = D - fs.nplanes < D, on the selected planes: every vertex of
# a leaf emits the keys of its d-faces (D-d of its planes, the face phase); the vertices of a
# face are those emitting its key, sorted.
function _add_faces!(cf::CutCellFaces{D},fs::CutFaces,planes) where D
  c, k = cf.cell, fs.nplanes
  index, pairs = cf.face_index, cf.face_vertex_pairs
  empty!(index)
  empty!(pairs)
  for n in get_leaves(c)
    phase = get_node_phase(c,n)
    for v in get_node_vertices(c,n)
      vplanes = get_vertex_planes(c,v)
      for mask in 1:(1<<D)-1
        count_ones(mask) == k || continue
        all(i -> !_inmask(mask,i) || _selected(planes,vplanes[i]),1:D) || continue
        key = (_key_planes(vplanes,mask),_face_phase(c,phase,vplanes,mask))
        f = get!(index,key) do
          append!(fs.face_planes,Iterators.take(key[1],k))
          push!(fs.face_phase,key[2])
          Int32(length(fs.face_phase))
        end
        push!(pairs,(f,v))
      end
    end
  end
  sort!(pairs;alg=QuickSort) # in place (the default sort allocates)
  for (i,(f,v)) in enumerate(pairs)
    i > 1 && pairs[i-1] == (f,v) && continue # a vertex of a face shared by two leaves
    push!(fs.face_vertices,v)
    f == length(fs.face_ptrs) || push!(fs.face_ptrs,length(fs.face_vertices))
  end
  isempty(pairs) || push!(fs.face_ptrs,length(fs.face_vertices)+1)
  return fs
end

_inmask(mask,i) = (mask >> (i-1)) & 1 == 1

# The planes of `vplanes` selected by `mask` (sorted), padded with zeros
function _key_planes(vplanes::NTuple{D,Int32},mask) where D
  ntuple(Val(D)) do j
    for i in 1:D
      _inmask(mask,i) && (j -= 1) == 0 && return vplanes[i]
    end
    zero(Int32)
  end
end

function _face_phase(c::CutCell{D,N},phase::NTuple{N,Int8},planes,mask) where {D,N}
  for i in 1:D
    p = planes[i]
    (_inmask(mask,i) && is_levelset_plane(c,p)) || continue
    phase = Base.setindex(phase,Int8(INTERFACE),get_plane_levelset(c,p))
  end
  phase
end

function _triangulate_face!(cf::CutCellFaces{D},fs::CutFaces,f::Integer) where D
  verts = cf.vertex_buffers[1]
  empty!(verts)
  append!(verts,view(fs.face_vertices,fs.face_ptrs[f]:fs.face_ptrs[f+1]-1))
  empty!(cf.plane_stack)
  k = fs.nplanes
  append!(cf.plane_stack,view(fs.face_planes,(f-1)*k+1:f*k))
  empty!(cf.apex_stack)
  _pull!(cf,fs,1,fs.d)
  push!(fs.face_simplex_ptrs,length(fs.simplex_vertices) ÷ (fs.d+1) + 1)
end

# Pulling triangulation of the dim-face with vertices `vertex_buffers[depth]` and planes
# `plane_stack`: cone from its smallest vertex over its facets not containing it. The facets
# are the planes of its vertices that are not its own; every such plane is a facet (simple
# polytopes), and the first vertex on it is the one that reports it.
function _pull!(cf::CutCellFaces,fs::CutFaces,depth::Integer,dim::Integer)
  c = cf.cell
  verts = cf.vertex_buffers[depth]
  if dim == 0
    @check length(verts) == 1
    append!(fs.simplex_vertices,cf.apex_stack)
    push!(fs.simplex_vertices,verts[1])
    return
  end
  v0 = minimum(verts)
  push!(cf.apex_stack,v0)
  v0planes = get_vertex_planes(c,v0)
  for (i,u) in enumerate(verts), p in get_vertex_planes(c,u)
    (p in cf.plane_stack || p in v0planes) && continue
    _reported_before(c,verts,i,p) && continue
    sub = cf.vertex_buffers[depth+1]
    empty!(sub)
    for w in verts
      p in get_vertex_planes(c,w) && push!(sub,w)
    end
    push!(cf.plane_stack,p)
    _pull!(cf,fs,depth+1,dim-1)
    pop!(cf.plane_stack)
  end
  pop!(cf.apex_stack)
  return
end

function _reported_before(c::CutCell,verts,i,p)
  for j in 1:i-1
    p in get_vertex_planes(c,verts[j]) && return true
  end
  false
end

## Display

function Base.show(io::IO,cf::CutCellFaces{D}) where D
  print(io,"CutCellFaces{$D}(faces per dimension: ",map(d -> num_faces(cf,d),0:D),")")
end
