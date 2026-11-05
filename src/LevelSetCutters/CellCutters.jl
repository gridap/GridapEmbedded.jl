# Copied from background-gmg, src/cutters/cellcutter/CellCutters.jl (the recursive cell cutter, originally
# PolytopalCutters.jl). Changes are marked with [changed].

"""
    CellCutter(p::Polytope,N::Integer;atol=1e-12,phases=ntuple(_ -> UNSET,N))

`Map` cutting the reference cell `p` by `N` affine level sets. Each level set is given by
its values at the vertices of `p`, e.g. the cell dof values of a P1 FE function, and the
result is a `CutCell`:

    cutter = CellCutter(TRI,2)
    cell_cuts = lazy_map(cutter,get_cell_dof_values(φ₁),get_cell_dof_values(φ₂))

Level sets are applied in order, so node phases are resolved from `φ₁` to `φN`. A level set
is considered to vanish on a vertex if its absolute value is at most `atol` (absolute: level
sets are meant to be distance functions) [changed: was relative to the largest value on
the cell].

Degenerate cuts (`D+1` planes through a point) throw a `DegenerateCutError`, unless the level
sets are shifted: with a first argument `shifts` (one `SVector{M}` per level set), level set
`k` is `φₖ + Σⱼ shifts[k][j] εʲ`, `ε → 0⁺` (see `CutCell`):

    cell_cuts = lazy_map(cutter,cell_shifts,φ₁_vals,φ₂_vals)

The zero tolerance of the ε-coefficients is also `atol` [changed: was relative].

`phases` selects the part of the tree that is built, one entry per level set: `IN` or `OUT`
(only that side of it is kept) or `UNSET = 0` (both sides; the default, all `UNSET`, builds
the whole tree). With every entry `IN` or `OUT` the tree is a single path. Nodes on an
excluded side are pruned, so the leaves are the selected phases only.
"""
struct CellCutter{D,N,R<:ReferenceCell{D}} <: Map
  refcell::R
  atol::Float64
  phases::NTuple{N,Int8}
end

function CellCutter(p::Polytope{D},N::Integer;atol::Real=1e-12,phases::Tuple=ntuple(_ -> UNSET,N)) where D
  @check length(phases) == N && all(∈((IN,OUT,UNSET)),phases) "phases: IN, OUT or UNSET per level set"
  refcell = ReferenceCell(p)
  CellCutter{D,N,typeof(refcell)}(refcell,atol,NTuple{N,Int8}(phases))
end

const Shifts{M} = AbstractVector{<:SVector{M}}

function Arrays.return_cache(k::CellCutter{D,N},shifts::Shifts{M},φs::Vararg{AbstractVector,N}) where {D,N,M}
  T = float(promote_type(map(eltype,φs)...))
  return CutCell{N,M+1}(k.refcell,T,k.atol)
end

function Arrays.evaluate!(c::CutCell,k::CellCutter{D,N},shifts::Shifts,φs::Vararg{AbstractVector,N}) where {D,N}
  @check c.refcell === k.refcell
  @check length(shifts) == N
  _init!(c,shifts,φs)
  _cut!(c,k.phases,ROOT,1)
  return c
end

# Unshifted level sets
_no_shifts(φs) = Fill(SVector{0,Float64}(),length(φs))
Arrays.return_cache(k::CellCutter{D,N},φs::Vararg{AbstractVector,N}) where {D,N} = return_cache(k,_no_shifts(φs),φs...)
Arrays.evaluate!(c::CutCell,k::CellCutter{D,N},φs::Vararg{AbstractVector,N}) where {D,N} = evaluate!(c,k,_no_shifts(φs),φs...)

# The root node is the cell: its vertices lie on their facets and take the given values.
function _init!(c::CutCell{D,N,T},shifts::Shifts,φs::NTuple{N,AbstractVector}) where {D,N,T}
  _reset!(c)
  r = c.refcell
  for k in 1:N
    @check length(φs[k]) == num_vertices(r)
    p = levelset_plane(c,k)
    c.planes[p] = levelset_plane(r,φs[k])
    c.offsets[p] = _epoly(T,c.planes[p].offset,shifts[k])
    c.tols[k] = _epoly(T,c.atol,map(_ -> c.atol,shifts[k])) # [changed: absolute tolerance]
  end
  for v in 1:num_vertices(r)
    x = VectorValue{D,T}(Tuple(r.vertex_coordinates[v]))
    values = ntuple(k -> T(φs[k][v]),Val(N))
    epoly = k -> _epoly(T,φs[k][v],shifts[k])
    _push_node_vertex!(c,_add_vertex!(c,r.vertex_to_facets[v],x,values,epoly))
  end
  _commit_node!(c,ntuple(_ -> Int8(UNSET),Val(N)))
  return c
end

# Resolve the phase of node `n` from level set `k` onwards. The node is cut by the first level
# set it is not on one side of, and the recursion continues on its (selected) children. A node
# on an excluded side of a level set is pruned.
function _cut!(c::CutCell{D,N},phases::NTuple{N,Int8},n::Int32,k::Integer) where {D,N}
  for l in k:N
    side = _side(c,n,l)
    if side == UNSET
      nin, nout = _split!(c,n,l,phases[l])
      iszero(nin) || _cut!(c,phases,nin,l+1)
      iszero(nout) || _cut!(c,phases,nout,l+1)
      return
    end
    _set_phase!(c,n,l,side)
    (phases[l] == UNSET || phases[l] == side) || return _prune!(c,n,l)
  end
end

_is_in(c::CutCell,v::Integer,k::Integer) = get_vertex_signs(c,v)[k] < 0

function _side(c::CutCell,n::Integer,k::Integer)
  has_in, has_out = false, false
  for v in get_node_vertices(c,n)
    _is_in(c,v,k) ? (has_in = true) : (has_out = true)
  end
  return (has_in && has_out) ? Int8(UNSET) : (has_in ? Int8(IN) : Int8(OUT))
end

# The selected children (`side`: IN, OUT or UNSET for both; 0 for a child not built) get the
# vertices on their side plus the vertices on the cut edges. With both, the latter are pushed
# to the IN child first, and copied over to the OUT child.
function _split!(c::CutCell,n::Int32,k::Integer,side::Int8)
  verts = get_node_vertices(c,n)
  phase = get_node_phase(c,n)
  nin = nout = zero(Int32)
  cut_vertices = 1:0
  if side != OUT
    for v in verts
      _is_in(c,v,k) && _push_node_vertex!(c,v)
    end
    cut_vertices = _push_cut_vertices!(c,verts,k)
    nin = _commit_node!(c,Base.setindex(phase,Int8(IN),k))
  end
  if side != IN
    for v in verts
      _is_in(c,v,k) || _push_node_vertex!(c,v)
    end
    if side == OUT
      _push_cut_vertices!(c,verts,k)
    else
      for i in cut_vertices
        _push_node_vertex!(c,c.node_vertices[i])
      end
    end
    nout = _commit_node!(c,Base.setindex(phase,Int8(OUT),k))
  end
  _set_children!(c,n,k,(nin,nout))
  return nin, nout
end

# Pushes the vertices where level set `k` cuts the edges of the node with vertices `verts`;
# returns their range in `node_vertices`. They lie on level set `k`, so they can only exist
# already if another node has been split by it (another branch of the tree).
function _push_cut_vertices!(c::CutCell,verts,k::Integer)
  start = length(c.node_vertices) + 1
  shared = c.levelset_split[k]
  for i in eachindex(verts), j in i+1:lastindex(verts)
    u, v = verts[i], verts[j]
    (_is_in(c,u,k) != _is_in(c,v,k)) || continue
    _are_adjacent(get_vertex_planes(c,u),get_vertex_planes(c,v)) || continue
    _push_node_vertex!(c,_cut_edge!(c,u,v,k,shared))
  end
  c.levelset_split[k] = true
  start:length(c.node_vertices)
end

# Vertex where level set `k` cuts the edge `(u,v)`. It lies on the planes of the edge and on
# the level set, and is created only if no other branch has created it already.
function _cut_edge!(c::CutCell{D,N},u::Integer,v::Integer,k::Integer,shared::Bool) where {D,N}
  planes = _edge_cut_planes(get_vertex_planes(c,u),get_vertex_planes(c,v),levelset_plane(c,k))
  if shared
    w = findfirst(isequal(planes),c.vertex_planes)
    isnothing(w) || return Int32(w)
  end
  x = compute_vertex_coordinates(c,planes)
  values = ntuple(l -> get_plane(c,levelset_plane(c,l))(x),Val(N))
  return _add_vertex!(c,planes,x,values,l -> _vertex_epoly(c,planes,l))
end

# The ε-polynomial of level set `l` on the vertex on `planes`. Only offsets are shifted, so
# the vertex is `x(ε) = Σⱼ xⱼ εʲ⁻¹` with `A xⱼ = -bⱼ` (`A` the normals of the planes, `bⱼ`
# the ε-coefficients of their offsets). It depends on the planes only, not on the branch that
# creates the vertex.
function _vertex_epoly(c::CutCell{D,N,T,R,K},planes::NTuple{D,Int32},l::Integer) where {D,N,T,R,K}
  B = SMatrix{D,K,T}(ntuple(i -> -c.offsets[planes[(i-1)%D+1]][(i-1)÷D+1],Val(D*K)))
  X = _normal_matrix(c,planes)\B
  p = levelset_plane(c,l)
  transpose(X)*SVector(Tuple(get_plane(c,p).normal)) + c.offsets[p]
end

## Plane-set combinatorics
# Vertices of a simple polytope share an edge iff they share D-1 planes.

_are_adjacent(a::NTuple{D},b::NTuple{D}) where D = count(p -> p ∈ b,a) == D-1

# Planes of the edge (a ∩ b) plus the cutting plane `p`. Level sets are applied in increasing
# order and come after the facets, so `p` is larger than all other planes and the result
# stays sorted.
function _edge_cut_planes(a::NTuple{D,Int32},b::NTuple{D,Int32},p::Int32) where D
  j = findfirst(x -> x ∉ b,a)::Int
  @check p > last(a) && p > last(b)
  return ntuple(i -> (i < j) ? a[i] : ((i < D) ? a[i+1] : p),Val(D))
end
