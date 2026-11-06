# Tests of the recursive cell cutter (src/LevelSetCutters/{CellCutters,CutFaces}.jl); helpers shared
# by its test modules. Geometric helpers, only used to check the topology against analytic
# results.

using GridapEmbedded.LevelSetCutters: Plane, ReferenceCell, CutCell, CellCutter,
  num_levelsets, num_planes, get_plane, get_plane_normal, levelset_plane, is_levelset_plane,
  get_plane_levelset, get_vertex_planes, compute_vertex_coordinates, get_vertex_values, get_vertex_signs, get_node_vertices,
  get_node_phase, get_node_children, get_node_levelset, get_node_vertices_on_plane, isleaf,
  ispruned,
  get_leaves, CutCellFaces, CellFaceExtractor, get_cut_cell, get_face_planes, get_face_phase,
  get_face_normal, num_simplices, get_simplex_vertices, get_face_simplices, get_face_quadrature
using GridapEmbedded.Interfaces: IN, OUT, INTERFACE, UNSET

using LinearAlgebra
using Random
using StaticArrays
using Gridap.TensorValues

# Volume of a convex node, as a sum of pyramids from its centroid to its facets.
function node_volume(c::CutCell{D},n) where D
  verts = get_node_vertices(c,n)
  x0 = sum(v -> get_vertex_coordinates(c,v),verts)/length(verts)
  planes = unique(p for v in verts for p in get_vertex_planes(c,v))
  vol = 0.0
  for p in planes
    xs = [get_vertex_coordinates(c,v) for v in get_node_vertices_on_plane(c,n,p)]
    @assert length(xs) ≥ D
    h = abs(get_plane_normal(c,p)⋅(first(xs)-x0))
    vol += h*facet_measure(xs,get_plane_normal(c,p))/D
  end
  return vol
end

facet_measure(xs::Vector{<:VectorValue{2}},n) = (@assert length(xs) == 2; norm(xs[2]-xs[1]))

function facet_measure(xs::Vector{<:VectorValue{3}},n)
  xc = sum(xs)/length(xs)
  r, i = findmax(x -> norm(x-xc),xs)
  iszero(r) && return 0.0 # all vertices coincide (shifted degeneracies)
  e1 = (xs[i]-xc)/r
  e2 = n × e1
  θs = map(x -> atan((x-xc)⋅e2,(x-xc)⋅e1),xs)
  ys = xs[sortperm(θs)]
  return sum(i -> norm((ys[i]-xc) × (ys[mod1(i+1,length(ys))]-xc))/2,eachindex(ys))
end

cell_volume(::CutCell{2}) = 1/2
cell_volume(::CutCell{3}) = 1/6

# Leaf with the given phase, or nothing if that region is empty.
function find_leaf(c,phase)
  i = findfirst(n -> get_node_phase(c,n) == Int8.(phase),collect(get_leaves(c)))
  return isnothing(i) ? nothing : collect(get_leaves(c))[i]
end

region_volume(c,phase) = (n = find_leaf(c,phase); isnothing(n) ? 0.0 : node_volume(c,n))

function cut(p,φs...;phases=map(_ -> UNSET,φs))
  k = CellCutter(p,length(φs);phases)
  c = return_cache(k,φs...)
  return evaluate!(c,k,φs...)
end

function cut(p,shifts::AbstractVector{<:SVector},φs...;phases=map(_ -> UNSET,φs))
  k = CellCutter(p,length(φs);phases)
  c = return_cache(k,shifts,φs...)
  return evaluate!(c,k,shifts,φs...)
end

# Vertex values of the affine function x -> a⋅x + b on the reference simplex
vertex_values(p,a,b) = map(x -> a⋅x + b,get_vertex_coordinates(p))

# Shifts of N level sets resolving every degeneracy but coincident planes: φₖ + εᵏ
independent_shifts(N) = [SVector(ntuple(j -> Float64(j == k),N)) for k in 1:N]

# Shifts translating the level set x -> a⋅x + b by the direction (ε,ε²,…). They resolve the
# degeneracies between the level sets and the cell, as long as no D+1 level sets meet at a
# point (e.g. the facets of a convex cell).
translation_shifts(a::VectorValue) = -SVector(Tuple(a))

_e(D,i) = VectorValue(ntuple(j -> Float64(i == j),D))

# N distinct planes a⋅x + b with a ∈ {eᵢ, eᵢ-eⱼ, 1} (up to sign) through points of the lattice
# {0,1/2,1}ᴰ: they meet the vertices, edges and facets of the cell, and each other there.
function degenerate_planes(rng,p,N)
  D = num_dims(p)
  normals = [[_e(D,i) for i in 1:D]; [_e(D,i)-_e(D,j) for i in 1:D for j in i+1:D]; [sum(i -> _e(D,i),1:D)]]
  planes = Set{Tuple{VectorValue{D,Float64},Float64}}()
  while length(planes) < N
    a = rand(rng,normals)
    q = VectorValue(ntuple(_ -> rand(rng,(0.0,0.5,1.0)),D))
    push!(planes,(a,-a⋅q))
  end
  map(((a,b),) -> (rand(rng,Bool) ? (a,b) : (-a,-b)),collect(planes))
end

# The facets a⋅x + b ≤ 0 of an axis-aligned box or a Kuhn simplex with vertices on the
# lattice (1/2)ℤᴰ, near the cell: background cells whose planes meet the vertices, edges and
# facets of the cell.
function lattice_cell_planes(rng,p)
  D = num_dims(p)
  e(i) = _e(D,i)
  q = VectorValue(ntuple(_ -> rand(rng,(-0.5,0.0,0.5)),D))
  h = rand(rng,(0.5,1.0))
  rand(rng,Bool) && return [[(-e(i),q[i]) for i in 1:D]; [(e(i),-q[i]-h) for i in 1:D]]
  # 0 ≤ y_σ₁ ≤ … ≤ y_σD ≤ h, y = x - q
  σ = randperm(rng,D)
  [[(-e(σ[1]),q[σ[1]])]; [(e(σ[i])-e(σ[i+1]),q[σ[i+1]]-q[σ[i]]) for i in 1:D-1]; [(e(σ[D]),-q[σ[D]]-h)]]
end
