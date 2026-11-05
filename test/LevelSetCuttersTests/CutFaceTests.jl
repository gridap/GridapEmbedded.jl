# Copied from background-gmg, test/cutters/CutFaceTests.jl.
module CutFaceTests

using Test
using Random
using LinearAlgebra
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using StaticArrays

include("CellCutterHelpers.jl")

function extract(c::CutCell)
  cf = return_cache(CellFaceExtractor(),c)
  evaluate!(cf,CellFaceExtractor(),c)
end

measure(cf,d,f) = sum(get_face_quadrature(cf,d,f,1)[2])
integrate(g,cf,d,f,degree) = ((xs,ws) = get_face_quadrature(cf,d,f,degree); sum(w*g(x) for (x,w) in zip(xs,ws)))

# Independent measure of a d-face from its vertices (convex)
function face_measure(cf::CutCellFaces{D},d,f) where D
  xs = [get_vertex_coordinates(cf,v) for v in get_face_vertices(cf,d,f)]
  d == 0 && return 1.0
  d == 1 && return norm(xs[2]-xs[1])
  d == D && return node_volume(get_cut_cell(cf),_leaf(cf,f))
  facet_measure(xs,get_face_normal(cf,f))
end
_leaf(cf,f) = collect(get_leaves(get_cut_cell(cf)))[f]

# The (D-1)-faces bounding leaf (D-face) f, with the outward normal
function leaf_facets(cf::CutCellFaces{D,N},f) where {D,N}
  c = get_cut_cell(cf)
  phase = get_face_phase(cf,D,f)
  verts = get_face_vertices(cf,D,f)
  out = Tuple{Int,VectorValue{D,Float64}}[]
  for g in 1:num_faces(cf,D-1)
    all(v -> v in verts,get_face_vertices(cf,D-1,g)) || continue
    p = only(get_face_planes(cf,D-1,g))
    gphase = get_face_phase(cf,D-1,g)
    k = is_levelset_plane(c,p) ? get_plane_levelset(c,p) : 0
    all(l -> l == k || gphase[l] == phase[l],1:N) || continue
    n = get_face_normal(cf,g)
    (k > 0 && phase[k] == OUT) && (n = -n)
    push!(out,(g,n))
  end
  out
end

# The facets of the D-simplices lying on a plane, and the (D-1)-simplices of the faces, as
# sorted vertex sets: the triangulation of the faces is that of the leaves (conforming)
function boundary_simplices(cf::CutCellFaces{D}) where D
  c = get_cut_cell(cf)
  bnd = Set{Vector{Int32}}()
  for s in 1:num_simplices(cf,D)
    S = collect(get_simplex_vertices(cf,D,s))
    for i in 1:D+1
      t = sort(deleteat!(copy(S),i))
      common = intersect(map(v -> Set(get_vertex_planes(c,v)),t)...)
      isempty(common) || push!(bnd,t)
    end
  end
  faces = Set(sort(collect(get_simplex_vertices(cf,D-1,s))) for s in 1:num_simplices(cf,D-1))
  bnd, faces
end

function test_faces(cf::CutCellFaces{D,N}) where {D,N}
  c = get_cut_cell(cf)
  # Euler characteristic of the cell (a ball)
  @test sum(d -> (-1)^d*num_faces(cf,d),0:D) == 1
  for d in 0:D, f in 1:num_faces(cf,d)
    planes = get_face_planes(cf,d,f)
    phase = get_face_phase(cf,d,f)
    verts = get_face_vertices(cf,d,f)
    @test length(planes) == D-d && issorted(planes)
    @test length(verts) ≥ d+1
    # its vertices are on its planes and on its side of every level set
    for v in verts
      @test all(p -> p in get_vertex_planes(c,v),planes)
      for k in 1:N
        σ = get_vertex_signs(c,v)[k]
        if phase[k] == INTERFACE
          @test levelset_plane(c,k) in planes
        else
          @test (phase[k] == IN ? σ ≤ 0 : σ ≥ 0)
        end
      end
    end
    # its simplices use its vertices and cover it
    for s in get_face_simplices(cf,d,f)
      @test all(v -> v in verts,get_simplex_vertices(cf,d,s))
    end
    @test measure(cf,d,f) ≈ face_measure(cf,d,f) atol=1e-12
  end
  # faces are unique
  for d in 0:D
    @test allunique([(collect(get_face_planes(cf,d,f)),get_face_phase(cf,d,f)) for f in 1:num_faces(cf,d)])
  end
  # the leaves partition the cell; polynomials integrated exactly
  @test sum(f -> measure(cf,D,f),1:num_faces(cf,D)) ≈ cell_volume(c)
  g(x) = D == 2 ? x[1]^2*x[2] : x[1]*x[2]*x[3]
  exact = D == 2 ? 2*1/factorial(5) : 1/factorial(6)
  @test sum(f -> integrate(g,cf,D,f,3),1:num_faces(cf,D)) ≈ exact
  # divergence theorem on every leaf, with the face normals
  v(x) = D == 2 ? VectorValue(x[1]^2+x[2],x[1]*x[2]) : VectorValue(x[1]^2+x[2],x[1]*x[2],x[3]*x[1])
  divv(x) = D == 2 ? 3x[1] : 4x[1]
  for f in 1:num_faces(cf,D)
    lhs = integrate(divv,cf,D,f,2)
    rhs = sum(((h,n),) -> integrate(x -> v(x)⋅n,cf,D-1,h,2),leaf_facets(cf,f))
    @test lhs ≈ rhs atol=1e-12
  end
  # conforming: the leaf triangulations restrict to the face triangulations
  bnd, faces = boundary_simplices(cf)
  @test bnd == faces
end

@testset "Faces of one cut cell" begin
  # x < 0.6, y < 0.7 on the triangle
  c = cut(TRI,vertex_values(TRI,VectorValue(1.0,0.0),-0.6),vertex_values(TRI,VectorValue(0.0,1.0),-0.7))
  cf = extract(c)
  @test map(d -> num_faces(cf,d),0:2) == [7,9,3]
  test_faces(cf)
  # interface pieces: φ₁ = 0 for y < 0.7 (length 0.4), φ₂ = 0 for x < 0.6 (length 0.3)
  Γ(k) = filter(f -> get_face_phase(cf,1,f)[k] == INTERFACE,1:num_faces(cf,1))
  @test sum(f -> measure(cf,1,f),Γ(1)) ≈ 0.4
  @test sum(f -> measure(cf,1,f),Γ(2)) ≈ 0.3
  @test all(f -> get_face_normal(cf,f) ≈ VectorValue(1.0,0.0),Γ(1))

  # x+y+z < 0.5 on the tetrahedron: interface triangle of area √3/2 0.5²
  c = cut(TET,vertex_values(TET,VectorValue(1.0,1.0,1.0),-0.5))
  cf = extract(c)
  test_faces(cf)
  Γ = filter(f -> get_face_phase(cf,2,f)[1] == INTERFACE,1:num_faces(cf,2))
  @test sum(f -> measure(cf,2,f),Γ) ≈ sqrt(3)/2*0.5^2

  # uncut cell: the faces of the simplex
  c = cut(TET,vertex_values(TET,VectorValue(1.0,0.0,0.0),-2.0))
  cf = extract(c)
  @test map(d -> num_faces(cf,d),0:3) == [4,6,4,1]
  test_faces(cf)
end

@testset "Faces of random cuts" begin
  rng = MersenneTwister(1)
  for p in (TRI,TET), N in 1:4, _ in 1:20
    D = num_dims(p)
    φs = map(1:N) do _
      a = VectorValue(ntuple(_ -> randn(rng),D))
      x = VectorValue(ntuple(_ -> rand(rng)/D,D))
      vertex_values(p,a,-a⋅x)
    end
    c = try cut(p,φs...) catch e
      e isa DegenerateCutError ? continue : rethrow()
    end
    test_faces(extract(c))
  end
end

@testset "Faces of random degenerate cuts" begin
  # Zero-measure faces included: the complex, its triangulation and the quadrature still hold
  rng = MersenneTwister(3)
  φs(p,planes) = [vertex_values(p,a,b) for (a,b) in planes]
  for p in (TRI,TET), N in 1:4, _ in 1:20
    planes = degenerate_planes(rng,p,N)
    test_faces(extract(cut(p,independent_shifts(N),φs(p,planes)...)))
  end
  for p in (TRI,TET), _ in 1:50
    planes = lattice_cell_planes(rng,p)
    test_faces(extract(cut(p,[translation_shifts(a) for (a,_) in planes],φs(p,planes)...)))
  end
end

# A face by its planes, phase, vertices and simplices (vertex ids of the cut cell)
_face_key(cf,d,f) = (collect(get_face_planes(cf,d,f)),get_face_phase(cf,d,f),sort(collect(get_face_vertices(cf,d,f))),
  sort([sort(collect(get_simplex_vertices(cf,d,s))) for s in get_face_simplices(cf,d,f)]))

@testset "Face selection" begin
  # The faces selected by dimension and planes are those of the full extraction
  rng = MersenneTwister(5)
  for p in (TRI,TET), N in 1:4, _ in 1:20
    D, F = num_dims(p), num_facets(p)
    planes = degenerate_planes(rng,p,N)
    c = cut(p,independent_shifts(N),[vertex_values(p,a,b) for (a,b) in planes]...)
    full = extract(c)
    for (dims,sel) in ((nothing,1:F),(0:D-1,nothing),([1],randsubseq(rng,1:F+N,0.5)),([D-1,D],[F+1]))
      k = CellFaceExtractor(;dims,planes=sel)
      cf = evaluate!(return_cache(k,c),k,c)
      for d in 0:D
        expected = (isnothing(dims) || d in dims) ?
          [_face_key(full,d,f) for f in 1:num_faces(full,d) if isnothing(sel) || all(∈(sel),get_face_planes(full,d,f))] : []
        @test sort([_face_key(cf,d,f) for f in 1:num_faces(cf,d)]) == sort(expected)
      end
    end
  end
end

@testset "Face extraction does not allocate" begin
  for (p,N) in ((TRI,2),(TET,3))
    D = num_dims(p)
    φs = [vertex_values(p,VectorValue(ntuple(j -> Float64(i == j),D)),-0.3) for i in 1:N]
    c = cut(p,φs...)
    cf = extract(c)
    k = CellFaceExtractor()
    evaluate!(cf,k,c)
    @test (@allocated evaluate!(cf,k,c)) == 0
  end
end

end # module
