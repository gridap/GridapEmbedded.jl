# Copied from background-gmg, test/cutters/CellCutterTests.jl.
module CellCutterTests

using Test
using Random
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.Geometry
using StaticArrays

include("CellCutterHelpers.jl")

@testset "ReferenceCell" begin
  for p in (TRI,TET)
    r = ReferenceCell(p)
    D = num_dims(p)
    for (v,x) in enumerate(r.vertex_coordinates)
      @test length(r.vertex_to_facets[v]) == D
      for f in 1:num_facets(r)
        g = r.facet_planes[f](x)
        (f ∈ r.vertex_to_facets[v]) ? (@test abs(g) < 1e-14) : (@test g < 0)
      end
      @test map(λ -> λ(x),r.barycentric) ≈ [i == v for i in 1:D+1]
    end
  end
end

@testset "One level set" begin
  # x < 0.3 on the triangle
  c = cut(TRI,vertex_values(TRI,VectorValue(1.0,0.0),-0.3))
  @test num_vertices(c) == 5
  @test num_nodes(c) == 3
  nin, nout = get_node_children(c,1)
  @test length(get_node_vertices(c,nin)) == 4
  @test length(get_node_vertices(c,nout)) == 3
  @test node_volume(c,nin) ≈ 0.3 - 0.3^2/2
  @test node_volume(c,nout) ≈ 0.7^2/2

  # x+y+z < 0.5 on the tetrahedron
  c = cut(TET,vertex_values(TET,VectorValue(1.0,1.0,1.0),-0.5))
  nin, nout = get_node_children(c,1)
  @test node_volume(c,nin) ≈ 0.5^3/6
  @test node_volume(c,nout) ≈ (1-0.5^3)/6

  # Uncut cell: single node, phase resolved on the root
  c = cut(TRI,vertex_values(TRI,VectorValue(1.0,0.0),-2.0))
  @test num_nodes(c) == 1
  @test get_node_phase(c,1) == (IN,)
  @test collect(get_leaves(c)) == [1]
end

@testset "Axis-aligned level sets" begin
  # x < a, y < b on the triangle
  a, b = 0.6, 0.7
  c = cut(TRI,vertex_values(TRI,VectorValue(1.0,0.0),-a),vertex_values(TRI,VectorValue(0.0,1.0),-b))
  @test region_volume(c,(IN,IN)) ≈ a*b - (a+b-1)^2/2
  @test region_volume(c,(IN,OUT)) ≈ (1-b)^2/2
  @test region_volume(c,(OUT,IN)) ≈ (1-a)^2/2
  @test isnothing(find_leaf(c,(OUT,OUT)))
  # The OUT side of φ₁ is entirely IN φ₂, so it is not split
  _, nout = get_node_children(c,1)
  @test isleaf(c,nout) && get_node_phase(c,nout) == (OUT,IN)

  # x < a, y < b, z < c on the tetrahedron, by inclusion-exclusion
  as = (0.5,0.4,0.45)
  φs = [vertex_values(TET,VectorValue(ntuple(j -> Float64(i == j),3)),-as[i]) for i in 1:3]
  c = cut(TET,φs...)
  tail(s) = max(s,0.0)^3/6
  vol = 1/6 - sum(a -> tail(1-a),as) +
    sum(tail(1-as[i]-as[j]) for i in 1:3 for j in i+1:3) - tail(1-sum(as))
  @test region_volume(c,(IN,IN,IN)) ≈ vol
  @test sum(n -> node_volume(c,n),get_leaves(c)) ≈ 1/6
end

# Invariants of the topology, for random level sets
function test_invariants(c::CutCell{D,N}) where {D,N}
  F = num_facets(c)
  # Vertices: unique keys, consistent with the planes and the level-set values
  @test allunique(map(v -> get_vertex_planes(c,v),1:num_vertices(c)))
  for v in 1:num_vertices(c)
    x = get_vertex_coordinates(c,v)
    @test x ≈ compute_vertex_coordinates(c,v) atol=1e-12 # stored = recomputed from the planes
    for p in get_vertex_planes(c,v)
      @test abs(get_plane(c,p)(x)) < 1e-10
    end
    for k in 1:N
      @test get_vertex_values(c,v)[k] ≈ get_plane(c,levelset_plane(c,k))(x) atol=1e-10
    end
    for f in 1:F
      @test get_plane(c,f)(x) < 1e-10
    end
  end
  # Signs: zero on the vertex's own level sets only, that of the value away from zero
  for v in 1:num_vertices(c), k in 1:N
    σ, φ = get_vertex_signs(c,v)[k], get_vertex_values(c,v)[k]
    @test iszero(σ) == (levelset_plane(c,k) in get_vertex_planes(c,v))
    abs(φ) > 1e-10 && @test σ == sign(φ)
  end
  # Tree
  for n in 1:num_nodes(c)
    phase = get_node_phase(c,n)
    for v in get_node_vertices(c,n), k in 1:N
      σ = get_vertex_signs(c,v)[k]
      (phase[k] == IN) && @test σ ≤ 0
      (phase[k] == OUT) && @test σ ≥ 0
    end
    isleaf(c,n) && continue
    if ispruned(c,n) # on an excluded side of the level set
      @test phase[-get_node_levelset(c,n)] != UNSET
      continue
    end
    k = get_node_levelset(c,n)
    @test phase[k] == UNSET
    nin, nout = get_node_children(c,n)
    @test !(iszero(nin) && iszero(nout))
    (iszero(nin) || iszero(nout)) && continue # one side selected
    @test node_volume(c,nin) + node_volume(c,nout) ≈ node_volume(c,n) atol=1e-12
    # Both children share the facet on the cutting plane
    p = levelset_plane(c,k)
    @test Set(get_node_vertices_on_plane(c,nin,p)) == Set(get_node_vertices_on_plane(c,nout,p))
  end
  # Leaves are regions of different phases; they partition the cell unless phases were selected
  leaves = collect(get_leaves(c))
  @test all(n -> UNSET ∉ get_node_phase(c,n),leaves)
  @test allunique(map(n -> get_node_phase(c,n),leaves))
  full = all(n -> isleaf(c,n) || (!ispruned(c,n) && all(!iszero,get_node_children(c,n))),1:num_nodes(c))
  full && @test sum(n -> node_volume(c,n),leaves) ≈ cell_volume(c)
end

@testset "Random level sets" begin
  rng = MersenneTwister(0)
  for p in (TRI,TET), N in 1:4, _ in 1:20
    D = num_dims(p)
    # Planes through random points of the cell, so that most of them cut it
    φs = map(1:N) do _
      a = VectorValue(ntuple(_ -> randn(rng),D))
      x = VectorValue(ntuple(_ -> rand(rng)/D,D))
      vertex_values(p,a,-a⋅x)
    end
    c = try cut(p,φs...) catch e
      e isa DegenerateCutError ? continue : rethrow()
    end
    test_invariants(c)
  end
end

@testset "Degeneracies" begin
  # Level set through a vertex of the cell
  @test_throws DegenerateCutError cut(TRI,[0.0,1.0,-1.0])
  # Two level sets meeting on a facet of the cell: φ₂ vanishes on Γ₁ ∩ {y = 0}
  φ1 = vertex_values(TRI,VectorValue(1.0,0.0),-0.3)
  φ2 = vertex_values(TRI,VectorValue(-1.0,1.0),0.3)
  @test_throws DegenerateCutError cut(TRI,φ1,φ2)
  # Three level sets meeting at a point in 2D
  φ3 = vertex_values(TRI,VectorValue(0.0,1.0),-0.2)
  φ4 = vertex_values(TRI,VectorValue(1.0,1.0),-0.5)
  @test_throws DegenerateCutError cut(TRI,φ1,φ3,φ4)
  # ... which is fine in 3D, where the triple point is a vertex
  ψs = [vertex_values(TET,VectorValue(ntuple(j -> Float64(i == j),3)),-0.2) for i in 1:3]
  c = cut(TET,ψs...)
  @test count(v -> all(is_levelset_plane(c,p) for p in get_vertex_planes(c,v)),1:num_vertices(c)) == 1
end

@testset "Shifted degeneracies" begin
  # Level set φ = x - y through a vertex: the vertex takes the side of the shift
  φ = [0.0,1.0,-1.0]
  for (s,σ) in ((1.0,OUT),(-1.0,IN))
    c = cut(TRI,[SVector(s)],φ)
    test_invariants(c)
    @test get_vertex_signs(c,1) == (Int8(σ == IN ? -1 : 1),)
    @test region_volume(c,(IN,)) ≈ 1/4
  end
  # ... unless the shift vanishes too
  @test_throws DegenerateCutError cut(TRI,[SVector(0.0)],φ)
  # Level set φ = y on the facet y = 0: no cut, or a zero-measure IN region
  φ = vertex_values(TRI,VectorValue(0.0,1.0),0.0)
  c = cut(TRI,[SVector(1.0)],φ)
  @test num_nodes(c) == 1 && get_node_phase(c,1) == (OUT,)
  c = cut(TRI,[SVector(-1.0)],φ)
  test_invariants(c)
  @test region_volume(c,(IN,)) ≈ 0 atol=1e-14
  # The degeneracies of the unshifted tests. Two level sets meeting on a facet: translating
  # them separates them from the facet
  a1, a2, a3, a4 = VectorValue(1.0,0.0), VectorValue(-1.0,1.0), VectorValue(0.0,1.0), VectorValue(1.0,1.0)
  φ1, φ2, φ3, φ4 = vertex_values(TRI,a1,-0.3), vertex_values(TRI,a2,0.3), vertex_values(TRI,a3,-0.2), vertex_values(TRI,a4,-0.5)
  test_invariants(cut(TRI,translation_shifts.([a1,a2]),φ1,φ2))
  # Three level sets meeting at a point: a common translation keeps them concurrent, shifting
  # them independently does not
  @test_throws DegenerateCutError cut(TRI,translation_shifts.([a1,a3,a4]),φ1,φ3,φ4)
  test_invariants(cut(TRI,independent_shifts(3),φ1,φ3,φ4))
end

# A neighbour sees the same planes with the opposite orientation (and shifts): opposite signs,
# the same leaves with opposite phases.
function test_opposite(p,shifts,φs)
  c, c2 = cut(p,shifts,φs...), cut(p,map(-,shifts),map(-,φs)...)
  test_invariants(c)
  for v in 1:num_dims(p)+1
    @test get_vertex_signs(c2,v) == map(-,get_vertex_signs(c,v))
  end
  for n in get_leaves(c)
    phase = map(s -> s == IN ? OUT : IN,get_node_phase(c,n))
    @test node_volume(c,n) ≈ region_volume(c2,phase) atol=1e-12
  end
end

@testset "Random degenerate level sets" begin
  rng = MersenneTwister(2)
  for p in (TRI,TET), N in 1:4, _ in 1:50
    planes = degenerate_planes(rng,p,N)
    test_opposite(p,independent_shifts(N),[vertex_values(p,a,b) for (a,b) in planes])
  end
end

@testset "Degenerate background cells" begin
  # The mesh intersection: translation shifts, level sets on lattice boxes and Kuhn simplices
  rng = MersenneTwister(3)
  for p in (TRI,TET), _ in 1:200
    planes = lattice_cell_planes(rng,p)
    test_opposite(p,[translation_shifts(a) for (a,_) in planes],[vertex_values(p,a,b) for (a,b) in planes])
  end
end

# A leaf by its phase and the plane keys of its vertices (vertex ids depend on the tree)
_leaf_key(c,n) = (get_node_phase(c,n),sort([get_vertex_planes(c,v) for v in get_node_vertices(c,n)]))

@testset "Phase selectors" begin
  # The leaves of a selected tree are the leaves of the full tree with the selected phases;
  # with every level set selected (IN or OUT) the tree is a single path
  rng = MersenneTwister(4)
  for p in (TRI,TET), N in 1:3, _ in 1:20
    planes = degenerate_planes(rng,p,N)
    shifts, φs = independent_shifts(N), [vertex_values(p,a,b) for (a,b) in planes]
    full = cut(p,shifts,φs...)
    for phases in Iterators.product(ntuple(_ -> (IN,OUT,UNSET),N)...)
      c = cut(p,shifts,φs...;phases)
      test_invariants(c)
      selected(phase) = all(k -> phases[k] == UNSET || phase[k] == phases[k],1:N)
      @test sort([_leaf_key(c,n) for n in get_leaves(c)]) ==
            sort([_leaf_key(full,n) for n in get_leaves(full) if selected(get_node_phase(full,n))])
      if UNSET ∉ phases
        @test all(n -> count(!iszero,get_node_children(c,n)) ≤ 1,1:num_nodes(c))
      end
    end
  end
  @test_throws Exception CellCutter(TRI,2;phases=(IN,))
end

@testset "Display" begin
  c = cut(TRI,vertex_values(TRI,VectorValue(1.0,0.0),-0.6),vertex_values(TRI,VectorValue(0.0,1.0),-0.7))
  s = sprint(show,MIME"text/plain"(),c)
  @test occursin("5 nodes, 3 leaves",s)
  @test occursin("(IN,OUT)",s)
  # The OUT side of φ₁ is IN φ₂: pruned when selecting (OUT,OUT)
  c = cut(TRI,vertex_values(TRI,VectorValue(1.0,0.0),-0.6),vertex_values(TRI,VectorValue(0.0,1.0),-0.7);phases=(OUT,OUT))
  s = sprint(show,MIME"text/plain"(),c)
  @test occursin("2 nodes, 0 leaves",s)
  @test occursin("pruned by φ2",s)
end

## Copied from background-gmg, test/cutters/MapTests.jl


function physical_volume(c::CutCell{D},n,cell_coords) where D
  x0 = cell_coords[1]
  J = [cell_coords[j+1][i] - x0[i] for i in 1:D, j in 1:D]
  return node_volume(c,n)*abs(det(J))
end

@testset "lazy_map over a mesh" begin
  model = simplexify(CartesianDiscreteModel((0,1,0,1),(4,4)))
  V = FESpace(model,ReferenceFE(lagrangian,Float64,1))
  # Affine level sets, so their P1 interpolants are exact
  φ1 = interpolate(x -> x[1] + 0.5*x[2] - 0.6,V)
  φ2 = interpolate(x -> x[2] - 0.35,V)

  cutter = CellCutter(TRI,2)
  cell_cuts = lazy_map(cutter,get_cell_dof_values(φ1),get_cell_dof_values(φ2))
  cell_coords = get_cell_coordinates(get_triangulation(model))

  vols = Dict(phase => 0.0 for phase in ((IN,IN),(IN,OUT),(OUT,IN),(OUT,OUT)))
  cache = array_cache(cell_cuts)
  for cell in eachindex(cell_cuts)
    c = getindex!(cache,cell_cuts,cell)
    for n in get_leaves(c)
      vols[get_node_phase(c,n)] += physical_volume(c,n,cell_coords[cell])
    end
  end
  # {x + y/2 < 0.6} ∩ {y < 0.35} in the unit square
  @test vols[(IN,IN)] ≈ 0.6*0.35 - 0.35^2/4
  @test vols[(IN,OUT)] ≈ 0.6*0.65 - (1 - 0.35^2)/4
  @test sum(values(vols)) ≈ 1
end

@testset "No allocations" begin
  for (p,N) in ((TRI,2),(TET,3))
    D = num_dims(p)
    φs = [vertex_values(p,VectorValue(ntuple(j -> Float64(i == j),D)),-0.3) for i in 1:N]
    k = CellCutter(p,N)
    c = return_cache(k,φs...)
    evaluate!(c,k,φs...)
    @test (@allocated evaluate!(c,k,φs...)) == 0
  end
end

end # module
