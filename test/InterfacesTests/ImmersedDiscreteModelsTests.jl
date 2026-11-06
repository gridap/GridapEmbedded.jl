module ImmersedDiscreteModelsTests

# A simplicial (immersed) model cut by the cells of a background model (Cartesian or
# simplicial): `cut(bgmodel, imodel) -> ImmersedDiscreteModel`, geometry only, on
# hand-made fitted meshes.

using Test
using LinearAlgebra
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using GridapEmbedded
using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: get_imodel, get_bgmodel, get_cutfaces, get_iface_to_cutfaces,
  get_cutface_to_iface, get_bgcell_to_cutfaces, get_bgcell_to_ifaces, select_cutfaces,
  SubFaceTriangulation

## Fitted meshes (simplicial, inside the unit box)

rot(θ) = x -> VectorValue(cos(θ)*x[1]-sin(θ)*x[2], sin(θ)*x[1]+cos(θ)*x[2])
rot(θ,ϕ) = x -> (y = rot(θ)(VectorValue(x[1],x[2])); z = rot(ϕ)(VectorValue(y[2],x[3])); VectorValue(y[1],z[1],z[2]))

# Square/cube [-1,1]^D onto the disk/ball (elliptical grid mapping), then scaled, rotated, shifted
square_to_disk(x::VectorValue{2}) = VectorValue(x[1]*sqrt(1-x[2]^2/2), x[2]*sqrt(1-x[1]^2/2))
square_to_disk(x::VectorValue{3}) = VectorValue(
  x[1]*sqrt(1-x[2]^2/2-x[3]^2/2+x[2]^2*x[3]^2/3),
  x[2]*sqrt(1-x[3]^2/2-x[1]^2/2+x[3]^2*x[1]^2/3),
  x[3]*sqrt(1-x[1]^2/2-x[2]^2/2+x[1]^2*x[2]^2/3))

function ball_model(D,n;r=0.31,θ=0.3)
  c = D == 2 ? Point(0.512,0.487) : Point(0.512,0.487,0.505)
  R = D == 2 ? rot(θ) : rot(θ,0.7θ)
  simplexify(CartesianDiscreteModel(Tuple(repeat([-1,1],D)),ntuple(_ -> n,D);map=x -> c + r*R(square_to_disk(x))))
end

function box_model(D,n;θ=0.4)
  c = D == 2 ? Point(0.23,0.19) : Point(0.23,0.19,0.21)
  R = D == 2 ? rot(θ) : rot(θ,0.5θ)
  simplexify(CartesianDiscreteModel(Tuple(repeat([0,0.5],D)),ntuple(_ -> n,D);map=x -> c + R(x)))
end

# Grid-aligned: every fitted vertex on background planes (the background has 8^D cells)
aligned_model(D) = simplexify(CartesianDiscreteModel(Tuple(repeat([0.25,0.75],D)),ntuple(_ -> 4,D)))

## Helpers

# d-dimensional measure of the simplex with vertices S (any ambient dimension), via QR
function simplex_measure(S)
  d = length(S)-1
  d == 0 && return 1.0
  D = length(first(S))
  J = [(S[j+1]-S[1])[i] for i in 1:D, j in 1:d]
  abs(det(qr(J).R))/factorial(d)
end

# Largest relative difference between the measure of every d-iface and those of its cutfaces
function iface_coverage_error(m,d)
  topo = get_grid_topology(get_imodel(m))
  X = get_vertex_coordinates(topo)
  f_to_v = get_faces(topo,d,0)
  cf = get_cutfaces(m,d)
  f_to_cf = get_iface_to_cutfaces(m,d)
  maximum(1:length(f_to_v)) do f
    mf = simplex_measure(X[f_to_v[f]])
    mc = sum(c -> simplex_measure(cf.point_to_coords[cf.cell_to_points[c]]),f_to_cf[f];init=0.0)
    abs(mf-mc)/mf
  end
end

# The reference coordinates in the background cell and in the iface give back the points
function glue_error(m,d)
  cf = get_cutfaces(m,d)
  bgmap = get_cell_map(get_grid(get_bgmodel(m)))
  topo = get_grid_topology(get_imodel(m))
  X = get_vertex_coordinates(topo)
  f_to_v = get_faces(topo,d,0)
  c_to_f = get_cutface_to_iface(m,d)
  maximum(1:length(cf.cell_to_bgcell)) do c
    Xf = X[f_to_v[c_to_f[c]]]
    maximum(cf.cell_to_points[c]) do p
      x = cf.point_to_coords[p]
      ξ = cf.point_to_ircoords[p]
      xi = Xf[1] + sum(k -> ξ[k]*(Xf[k+1]-Xf[1]),1:d;init=zero(x))
      max(norm(bgmap[cf.cell_to_bgcell[c]](cf.point_to_rcoords[p])-x),norm(xi-x))
    end
  end
end

## Tests

cases = [
  ("2D ball",    ball_model(2,6),   CartesianDiscreteModel((0,1,0,1),(9,9))),
  ("2D box",     box_model(2,5),    CartesianDiscreteModel((0,1,0,1),(7,7))),
  ("2D aligned", aligned_model(2),  CartesianDiscreteModel((0,1,0,1),(8,8))),
  ("3D ball",    ball_model(3,3),   CartesianDiscreteModel((0,1,0,1,0,1),(5,5,5))),
  ("3D aligned", aligned_model(3),  CartesianDiscreteModel((0,1,0,1,0,1),(8,8,8))),
]

for (name,fit,bg) in cases, (bgname,bgm) in (("Cartesian",bg),("simplicial",simplexify(bg)))
  @testset "$name, $bgname background" begin
    D = num_cell_dims(fit)
    m = cut(bgm,fit)
    @test m isa ImmersedDiscreteModel
    @test get_imodel(m) === fit && get_bgmodel(m) === bgm
    # behaves as the immersed model
    @test num_cells(m) == num_cells(fit)
    @test get_grid_topology(m) === get_grid_topology(fit)

    # every d-iface is covered exactly once by its cutfaces, the vertices once each
    @test length(get_cutfaces(m,0).cell_to_bgcell) == num_vertices(fit)
    @test get_cutfaces(m,0).point_to_coords[reduce(vcat,get_iface_to_cutfaces(m,0))] ≈ get_vertex_coordinates(get_grid_topology(fit))
    @test all(d -> iface_coverage_error(m,d) < 1e-11,1:D)
    @test all(d -> glue_error(m,d) < 1e-12,0:D)

    # integration on the cutfaces, glued to the background
    for d in 1:D
      @test sum(∫(1)Measure(SubFaceTriangulation(m,d),1)) ≈ sum(∫(1)Measure(Triangulation(ReferenceFE{d},fit),1)) rtol=1e-12
    end
    f(x) = 1 + x[1] - 2x[2]^2
    @test sum(∫(f)Measure(SubFaceTriangulation(m,D),2)) ≈ sum(∫(f)Measure(Triangulation(fit),2)) rtol=1e-12

    # accessors
    for d in 0:D
      c_to_f = get_cutface_to_iface(m,d)
      f_to_c = get_iface_to_cutfaces(m,d)
      @test all(f -> all(c -> c_to_f[c] == f,f_to_c[f]),1:length(f_to_c))
      b_to_c = get_bgcell_to_cutfaces(m,d)
      cf = get_cutfaces(m,d)
      @test sum(length,b_to_c) == length(cf.cell_to_bgcell)
      @test all(K -> all(c -> cf.cell_to_bgcell[c] == K,b_to_c[K]),1:num_cells(bgm))
      b_to_f = get_bgcell_to_ifaces(m,d)
      @test all(K -> b_to_f[K] == sort(unique(c_to_f[b_to_c[K]])),1:num_cells(bgm))
    end
    # selection and restriction (only the points of the selected cutfaces are kept)
    is_cut = falses(num_cells(bgm)); is_cut[1:2:end] .= true
    ids = select_cutfaces(m,D;bgcell_mask=is_cut)
    @test all(c -> is_cut[get_cutfaces(m,D).cell_to_bgcell[c]],ids)
    st = ImmersedFaceData(get_cutfaces(m,D),ids)
    @test sort(unique(st.cell_to_points.data)) == 1:length(st.point_to_coords) == 1:length(st.point_to_ircoords)
    @test st.point_to_coords[st.cell_to_points.data] == get_cutfaces(m,D).point_to_coords[get_cutfaces(m,D).cell_to_points[ids].data]

    # the cut of the background by ∂Ω as an EmbeddedDiscretization
    cutgeo = EmbeddedDiscretization(m)
    ioc = only(cutgeo.ls_to_bgcell_to_inoutcut)
    @test count(==(CUT),ioc) > 0 && count(==(IN),ioc) > 0 && count(==(OUT),ioc) > 0
    vol = sum(∫(1)Measure(Triangulation(fit),1))
    @test sum(∫(1)Measure(Triangulation(cutgeo,PHYSICAL),1)) ≈ vol rtol=1e-12
    @test sum(∫(1)Measure(Triangulation(cutgeo,ACTIVE),1)) ≥ vol
    Γ = EmbeddedBoundary(cutgeo)
    nΓ = get_normal_vector(Γ)
    dΓ = Measure(Γ,2)
    @test sum(∫(1)dΓ) ≈ sum(∫(1)Measure(BoundaryTriangulation(fit),1)) rtol=1e-12
    @test norm(sum(∫(nΓ)dΓ)) < 1e-12
    @test sum(∫(CellField(x -> x,Γ)⋅nΓ)dΓ) ≈ D*vol rtol=1e-12 # divergence theorem: outward normals
  end
end

# The default cutter
fit, bg = ball_model(2,4), CartesianDiscreteModel((0,1,0,1),(6,6))
m1, m2 = cut(bg,fit), cut(RecursivePlaneCutter(),bg,fit)
@test get_cutfaces(m1,2).cell_to_points == get_cutfaces(m2,2).cell_to_points
@test get_cutfaces(m1,2).point_to_coords == get_cutfaces(m2,2).point_to_coords

end # module
