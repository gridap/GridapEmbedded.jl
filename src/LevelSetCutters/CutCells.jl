# Copied from background-gmg, src/cutters/cellcutter/{Planes,ReferenceCells,CutCells,Geometry}.jl
# (the recursive cell cutter, originally PolytopalCutters.jl), merged into one file. Changes are
# marked with [changed].

## Affine planes

"""
    struct Plane{D,T}

Affine function `g(x) = n⋅x + c` in reference coordinates. The plane is its zero set, and the
normal `n` points towards `g > 0`.

Cell facets are oriented outwards (the cell is `g ≤ 0`) and level sets towards `φ > 0`
(`IN` is `φ < 0`).
"""
struct Plane{D,T}
  normal::VectorValue{D,T}
  offset::T
end

(g::Plane)(x::VectorValue) = g.normal ⋅ x + g.offset

Base.zero(::Type{Plane{D,T}}) where {D,T} = Plane(zero(VectorValue{D,T}),zero(T))

Base.convert(::Type{Plane{D,T}},g::Plane{D}) where {D,T} = Plane{D,T}(g.normal,g.offset)

# Σᵢ αᵢ gᵢ, as Gridap's `linear_combination(coefficients, fields)`
function Fields.linear_combination(αs::AbstractVector,gs::AbstractVector{<:Plane})
  @check length(gs) == length(αs)
  normal = sum(i -> αs[i]*gs[i].normal, eachindex(gs))
  offset = sum(i -> αs[i]*gs[i].offset, eachindex(gs))
  return Plane(normal,offset)
end

## Reference cells

"""
    struct ReferenceCell{D,P}

Everything the cutter needs from the reference polytope `P`. It is independent of the
level sets, so it is built once per cell type.

- `vertex_to_facets[v]`: sorted facets of vertex `v`, its plane ids in the `CutCell`.
- `facet_planes[f]`: facet `f` as an outward-oriented `Plane`.
- `barycentric[v]`: the barycentric coordinate `λᵥ` as a `Plane`. A level set with vertex
  values `φ` is the affine function `Σᵥ φᵥ λᵥ`.

Only simplices are supported: the level sets are assumed affine on the cell, which is only
true for P1 level sets on simplices.
"""
struct ReferenceCell{D,P<:Polytope{D}}
  polytope::P
  vertex_coordinates::Vector{VectorValue{D,Float64}}
  vertex_to_facets::Vector{NTuple{D,Int32}}
  facet_planes::Vector{Plane{D,Float64}}
  barycentric::Vector{Plane{D,Float64}}
end

function ReferenceCell(p::Polytope{D}) where D
  @notimplementedif !is_simplex(p) "Level sets are only affine on simplices"
  xs = collect(get_vertex_coordinates(p))
  vertex_to_facets = map(fs -> NTuple{D,Int32}(sort(fs)), get_faces(p,0,D-1))
  facet_planes = map(get_facet_normal(p),get_faces(p,D-1,0)) do n, vs
    Plane(n,-(n⋅xs[first(vs)]))
  end
  barycentric = _barycentric_planes(xs)
  return ReferenceCell(p,xs,vertex_to_facets,facet_planes,barycentric)
end

# λᵢ(x) = aᵢ⋅x + bᵢ with λᵢ(xⱼ) = δᵢⱼ, i.e. the columns of inv([xⱼᵀ 1]).
function _barycentric_planes(xs::Vector{VectorValue{D,T}}) where {D,T}
  V = [(j ≤ D) ? xs[i][j] : one(T) for i in 1:D+1, j in 1:D+1]
  C = inv(V)
  return [Plane(VectorValue(ntuple(d -> C[d,i],Val(D))),C[D+1,i]) for i in 1:D+1]
end

num_vertices(r::ReferenceCell) = length(r.vertex_coordinates)
num_facets(r::ReferenceCell) = length(r.facet_planes)

# Affine function with the given vertex values
levelset_plane(r::ReferenceCell,φ::AbstractVector) = Fields.linear_combination(φ,r.barycentric)

## Phases
# GridapEmbedded's `IN`, `OUT` (and `INTERFACE` for the faces on a level set), stored as
# `Int8`. Unknown phases are Gridap's `UNSET`.

"""
    struct CutCell{D,N,T}

Topology of a reference cell cut by `N` affine level sets. It is the cache and the output
of `CellCutter`, so it is overwritten for every cell.

# Planes

Planes `1:F` are the facets of the cell and planes `F+k` the level sets `k = 1:N`, see
`levelset_plane`. All of them are stored as `Plane`s, so the geometry can be recovered.

# Vertices

A vertex is identified by the sorted ids of the `D` planes it lies on (`get_vertex_planes`).
This key is unique: vertices created in different branches of the tree are merged. Each
vertex also carries its coordinates (`get_vertex_coordinates`), the values of the `N` level
sets on it (`get_vertex_values`) and their signs, computed from its planes when the vertex is
created. Values of the level sets it lies on are exactly zero. The coordinates are a function
of the planes only, so they can be recomputed from them at any time
(`compute_vertex_coordinates`), e.g. after the level sets change without changing the
topology.

# Perturbation

Level set `k` may be shifted symbolically, `φₖ + Σⱼ sₖⱼ εʲ` with `ε → 0⁺` (`sₖ` given to
`CellCutter`). Signs are those of the shifted level sets: a value within the tolerance of
zero takes the sign of the first ε-coefficient that is not; the ε-coefficients are only
computed for the values within the tolerance (none, for a non-degenerate cut). Only offsets
move, so every
vertex `x(ε)` and every value is a polynomial in `ε`, and the cut is that of the shifted,
non-degenerate arrangement. Geometry is that of `ε = 0`: coincident vertices with different
keys and zero-measure nodes are possible.

# Tree

Each node is a convex polytope, given by its vertices (`get_node_vertices`). Node `1` is the
cell. A node cut by level set `k` (`get_node_levelset`) has two children `(in,out)`
(`get_node_children`), on the `φₖ < 0` and `φₖ > 0` sides.

The phase of a node (`get_node_phase`) gives its side (`IN`/`OUT`) of each level set, and is
`UNSET` where it is not known at that node. It is fully set on the leaves (`get_leaves`),
which partition the cell.

With a phase selector (`CellCutter(...; phases)`), only the selected children are built (the
other one is `0`), and a node lying on an excluded side of level set `k` is pruned
(`get_node_levelset` is `-k`): the leaves are then the selected part of the cell.

The edges and faces of a node are not stored: since nodes are simple polytopes, two vertices
share an edge iff they share `D-1` planes, and a node's facet on a plane is the set of its
vertices on that plane (`get_node_vertices_on_plane`).
"""
struct CutCell{D,N,T,R<:ReferenceCell{D},K}
  refcell::R
  atol::Float64
  planes::Vector{Plane{D,T}}
  # ε-polynomials (K coefficients, from ε⁰): plane offsets, zero tolerances of the level sets
  offsets::Vector{SVector{K,T}}
  tols::Vector{SVector{K,T}}
  # Vertices
  vertex_planes::Vector{NTuple{D,Int32}}
  vertex_coordinates::Vector{VectorValue{D,T}}
  vertex_values::Vector{NTuple{N,T}}
  vertex_signs::Vector{NTuple{N,Int8}}
  # Tree
  node_ptrs::Vector{Int32}
  node_vertices::Vector{Int32}
  node_phase::Vector{NTuple{N,Int8}}
  node_children::Vector{NTuple{2,Int32}}
  node_levelset::Vector{Int8}
  # Level sets that have split a node already (only their vertices can be found again)
  levelset_split::Vector{Bool}
end

function CutCell{N,K}(refcell::ReferenceCell{D},::Type{T},atol::Real) where {D,N,K,T}
  F = num_facets(refcell)
  planes = Vector{Plane{D,T}}(undef,F+N)
  copyto!(planes,refcell.facet_planes)
  offsets = Vector{SVector{K,T}}(undef,F+N)
  copyto!(offsets,map(g -> _epoly(T,g.offset,zero(SVector{K-1,T})),refcell.facet_planes))
  CutCell{D,N,T,typeof(refcell),K}(
    refcell, atol, planes, offsets, Vector{SVector{K,T}}(undef,N),
    NTuple{D,Int32}[], VectorValue{D,T}[], NTuple{N,T}[], NTuple{N,Int8}[],
    Int32[1], Int32[], NTuple{N,Int8}[], NTuple{2,Int32}[], Int8[], fill(false,N)
  )
end

const ROOT = Int32(1)

## Planes

num_levelsets(::CutCell{D,N}) where {D,N} = N
num_facets(c::CutCell) = num_facets(c.refcell)
num_planes(c::CutCell) = length(c.planes)

get_plane(c::CutCell,p::Integer) = c.planes[p]

levelset_plane(c::CutCell,k::Integer) = Int32(num_facets(c) + k)
is_levelset_plane(c::CutCell,p::Integer) = p > num_facets(c)
get_plane_levelset(c::CutCell,p::Integer) = (@check is_levelset_plane(c,p); p - num_facets(c))

## Vertices

num_vertices(c::CutCell) = length(c.vertex_planes)

get_vertex_planes(c::CutCell,v::Integer) = c.vertex_planes[v]
get_vertex_coordinates(c::CutCell,v::Integer) = c.vertex_coordinates[v]
get_vertex_values(c::CutCell,v::Integer) = c.vertex_values[v]
get_vertex_signs(c::CutCell,v::Integer) = c.vertex_signs[v]

## Tree

num_nodes(c::CutCell) = length(c.node_phase)

get_node_vertices(c::CutCell,n::Integer) = view(c.node_vertices,c.node_ptrs[n]:c.node_ptrs[n+1]-1)
get_node_phase(c::CutCell,n::Integer) = c.node_phase[n]
get_node_children(c::CutCell,n::Integer) = c.node_children[n]
get_node_levelset(c::CutCell,n::Integer) = Int(c.node_levelset[n])

isleaf(c::CutCell,n::Integer) = iszero(c.node_levelset[n])
ispruned(c::CutCell,n::Integer) = c.node_levelset[n] < 0
get_leaves(c::CutCell) = Iterators.filter(n -> isleaf(c,n), 1:num_nodes(c))

function get_node_vertices_on_plane(c::CutCell,n::Integer,p::Integer)
  Iterators.filter(v -> p ∈ get_vertex_planes(c,v), get_node_vertices(c,n))
end

## Construction
# The tree is stored as a table: the vertices of a node are pushed one by one, then the node
# is committed. Nodes are only ever appended.

function _reset!(c::CutCell)
  empty!(c.vertex_planes)
  empty!(c.vertex_coordinates)
  empty!(c.vertex_values)
  empty!(c.vertex_signs)
  resize!(c.node_ptrs,1)
  empty!(c.node_vertices)
  empty!(c.node_phase)
  empty!(c.node_children)
  empty!(c.node_levelset)
  fill!(c.levelset_split,false)
  return c
end

# `x`, `values`: the coordinates and level-set values of the vertex (ε = 0); `epoly(k)`: the
# ε-polynomial of level set k on it, only evaluated when `values[k]` is within the tolerance.
# The level sets it lies on are zero.
function _add_vertex!(c::CutCell{D,N,T},planes::NTuple{D,Int32},x::VectorValue{D},values::NTuple{N},epoly::F) where {D,N,T,F}
  @check all(i -> planes[i] < planes[i+1],1:D-1) "Plane ids must be sorted and unique"
  on = ntuple(k -> levelset_plane(c,k) ∈ planes,Val(N))
  signs = ntuple(k -> on[k] ? Int8(0) : _sign(values[k],c.tols[k],epoly,k),Val(N))
  for k in 1:N
    on[k] || !iszero(signs[k]) || throw(DegenerateCutError(planes,k))
  end
  push!(c.vertex_planes,planes)
  push!(c.vertex_coordinates,x)
  push!(c.vertex_values,ntuple(k -> on[k] ? zero(T) : T(values[k]),Val(N)))
  push!(c.vertex_signs,signs)
  return Int32(num_vertices(c))
end

_push_node_vertex!(c::CutCell,v::Integer) = push!(c.node_vertices,v)

function _commit_node!(c::CutCell{D,N},phase::NTuple{N,Int8}) where {D,N}
  push!(c.node_ptrs,length(c.node_vertices)+1)
  push!(c.node_phase,phase)
  push!(c.node_children,(zero(Int32),zero(Int32)))
  push!(c.node_levelset,zero(Int8))
  return Int32(num_nodes(c))
end

function _set_phase!(c::CutCell,n::Integer,k::Integer,side::Int8)
  c.node_phase[n] = Base.setindex(c.node_phase[n],side,k)
end

function _set_children!(c::CutCell,n::Integer,k::Integer,children::NTuple{2,Int32})
  c.node_children[n] = children
  c.node_levelset[n] = k
end

_prune!(c::CutCell,n::Integer,k::Integer) = (c.node_levelset[n] = -k; nothing)

## ε-polynomials
# Coefficients of ε⁰, ε¹, … as an `SVector`.

_epoly(::Type{T},a,s::SVector{M}) where {T,M} = SVector{M+1,T}(a,s...)

# Sign of level set k for ε → 0⁺: that of the first coefficient of its ε-polynomial above the
# tolerance. `v` is the ε⁰ coefficient; the others (`epoly(k)`) only if `v` is within it.
function _sign(v,tol::SVector{K},epoly::F,k) where {K,F}
  abs(v) > tol[1] && return Int8(sign(v))
  K == 1 && return Int8(0)
  p = epoly(k)
  for j in 2:K
    abs(p[j]) > tol[j] && return Int8(sign(p[j]))
  end
  return Int8(0)
end

## Non-degeneracy

"""
    struct DegenerateCutError <: Exception

Thrown when level set `levelset` (almost) vanishes on the vertex lying on `planes`, i.e.
when `D+1` planes (almost) meet at a point, also after the perturbation (if any).
"""
struct DegenerateCutError <: Exception
  planes::Tuple
  levelset::Int
end

function Base.showerror(io::IO,e::DegenerateCutError)
  print(io,"DegenerateCutError: level set $(e.levelset) vanishes on the vertex on planes ",
    "$(Int.(e.planes)), so more than D planes meet at a point. Shift the level sets.")
end

## Display

function Base.show(io::IO,c::CutCell{D,N}) where {D,N}
  print(io,"CutCell{$D,$N}($(num_vertices(c)) vertices, $(num_nodes(c)) nodes)")
end

function Base.show(io::IO,::MIME"text/plain",c::CutCell{D,N}) where {D,N}
  nleaves = count(_ -> true, get_leaves(c))
  print(io,"CutCell of a $(c.refcell.polytope) by $N level set(s): ",
    "$(num_vertices(c)) vertices, $(num_nodes(c)) nodes, $nleaves leaves")
  num_nodes(c) > 0 && _show_node(io,c,ROOT,"")
end

_phase_string(phase) = join(map(s -> s == IN ? "IN" : (s == OUT ? "OUT" : "·"), phase), ",")

function _show_node(io::IO,c::CutCell,n,prefix)
  nv = length(get_node_vertices(c,n))
  print(io,"\n",prefix,n,": (",_phase_string(get_node_phase(c,n)),") with $nv vertices")
  isleaf(c,n) && return
  ispruned(c,n) && return print(io,", pruned by φ",-get_node_levelset(c,n))
  print(io,", cut by φ",get_node_levelset(c,n))
  child_prefix = replace(prefix,"├─ " => "│  ","└─ " => "   ")
  nin, nout = get_node_children(c,n)
  iszero(nin) || _show_node(io,c,nin,child_prefix*(iszero(nout) ? "└─ " : "├─ "))
  iszero(nout) || _show_node(io,c,nout,child_prefix*"└─ ")
end

## Geometry
#
# Geometry is recovered from the topology. All quantities live in reference coordinates.

"""
    compute_vertex_coordinates(c::CutCell,v::Integer)
    compute_vertex_coordinates(c::CutCell,planes::NTuple{D,Int32})

Reference coordinates of vertex `v` (or of the vertex on `planes`), computed as the
intersection of its `D` planes. They only depend on the planes, which in turn depend affinely
on the vertex values of the level sets, so they can be recomputed (or differentiated) for new
level-set values with the same topology. They are those of the unshifted level sets
(`ε = 0`). `get_vertex_coordinates` returns the coordinates stored when the vertex was created.
"""
compute_vertex_coordinates(c::CutCell,v::Integer) = compute_vertex_coordinates(c,get_vertex_planes(c,v))

function compute_vertex_coordinates(c::CutCell{D,N,T},planes::NTuple{D,Int32}) where {D,N,T}
  b = SVector{D,T}(map(p -> -get_plane(c,p).offset,planes))
  return VectorValue(Tuple(_normal_matrix(c,planes)\b))
end

# Rows: the normals of `planes`
function _normal_matrix(c::CutCell{D,N,T},planes::NTuple{D}) where {D,N,T}
  SMatrix{D,D,T}(ntuple(i -> get_plane(c,planes[(i-1)%D+1]).normal[(i-1)÷D+1],Val(D*D)))
end

"""
    get_plane_normal(c::CutCell,p::Integer)

Unit normal of plane `p` in reference coordinates. It points outwards for cell facets and
towards `φ > 0` for level sets.
"""
function get_plane_normal(c::CutCell,p::Integer)
  n = get_plane(c,p).normal
  return n/norm(n)
end
