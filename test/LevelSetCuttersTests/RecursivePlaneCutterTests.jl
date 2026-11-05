module RecursivePlaneCutterTests

# RecursivePlaneCutter against LevelSetCutter: both cut the same piecewise-linear level sets
# (simplexified cells, values at the vertices) exactly, so labels agree and measures agree to
# round-off, for every CSG combination; the pieces themselves differ.

using Test
using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using GridapEmbedded
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.CSG: get_tree
using AbstractTrees: Leaves

rotz(θ) = x -> VectorValue(cos(θ)*x[1]-sin(θ)*x[2], sin(θ)*x[1]+cos(θ)*x[2], x[3])
flower(x) = (r=sqrt((x[1]-0.48)^2+(x[2]-0.51)^2); θ=atan(x[2]-0.51,x[1]-0.48); r-0.3-0.07cos(5θ))
popc(x) = (y = 3.2*(x-Point(0.5,0.5,0.5)); y⋅y - 0.6 - 0.25*(y[1]^4+y[2]^4+y[3]^4) + 0.15*cos(4y[1])*cos(4y[2]))

# (name, model, geometry, closed interface[, atol])
# Not closed where level sets cross: both cutters leave round-off pieces (~1e-16) there.
# atol = 0 where values vanish at nodes up to round-off: the default atol makes them IN, while
# LevelSetCutter takes the round-off sign, so the labels would differ (measures do not).
function cases()
  quad = CartesianDiscreteModel((0,1,0,1),(20,20))
  quad_odd = CartesianDiscreteModel((-1,2,-1,2),(31,31))
  tri = simplexify(CartesianDiscreteModel((0,1,0,1),(16,16)))
  hex = CartesianDiscreteModel((0,1,0,1,0,1),(8,8,8))
  hex_rot = CartesianDiscreteModel((0,1,0,1,0,1),(7,7,7);map=rotz(0.3))
  tet = simplexify(CartesianDiscreteModel((0,1,0,1,0,1),(6,6,6)))
  d1 = disk(0.85,x0=Point(0,0),name="d1"); d2 = disk(0.85,x0=Point(1,1),name="d2"); d3 = disk(0.85,x0=Point(1,0),name="d3")
  s1 = sphere(0.33,x0=Point(0.47,0.52,0.49),name="s1")
  s2 = sphere(0.25,x0=Point(0.62,0.40,0.55),name="s2")
  dn = doughnut(0.25,0.09,x0=Point(0.5,0.52,0.47),name="doughnut")
  cb = cube(L=0.41,x0=Point(0.52,0.47,0.51),name="cube")
  V2 = FESpace(quad,ReferenceFE(lagrangian,Float64,1))
  V3 = FESpace(tet,ReferenceFE(lagrangian,Float64,1))
  [
    ("quad_disk",        quad,     disk(0.31,x0=Point(0.47,0.53),name="d"),                      true),
    ("quad_three_disks", quad_odd, union(union(d1,d2),d3),                                        false),
    ("tri_setdiff",      tri,      setdiff(disk(0.35,x0=Point(0.5,0.5),name="a"),disk(0.2,x0=Point(0.62,0.55),name="b")), true),
    ("tri_intersect",    tri,      intersect(disk(0.35,x0=Point(0.45,0.5),name="a"),disk(0.3,x0=Point(0.62,0.55),name="b")), true),
    ("quad_discrete",    quad,     DiscreteGeometry(interpolate(flower,V2),quad,name="flower"),   true),
    ("quad_nodes_line",  quad,     AnalyticalGeometry(x->x[1]-0.5,name="line"),                   false), # zero at nodes
    ("quad_nodes_diag",  quad,     AnalyticalGeometry(x->x[1]-x[2]-0.1,name="diag"),              false, 0.0),
    ("hex_sphere",       hex,      s1,                                                            true),
    ("hex_rot_multi",    hex_rot,  union(dn,!s2),                                                 false),
    ("hex_cube",         hex,      cb,                                                            true),  # 6 level sets
    ("tet_setdiff",      tet,      setdiff(s1,s2),                                                true),
    ("tet_discrete",     tet,      DiscreteGeometry(interpolate(popc,V3),tet,name="popcorn"),     true),
    ("hex_nodes_plane",  hex,      plane(x0=Point(0.5,0.5,0.5),v=VectorValue(1.,0.,0.),name="plane"), false),
  ]
end

leaf_names(geo) = unique([l.data[2] for l in Leaves(get_tree(geo))])

function boundary_quantities(Γ)
  dΓ = Measure(Γ,2)
  nΓ = get_normal_vector(Γ)
  (sum(∫(1)dΓ), sum(∫(nΓ)dΓ), sum(∫(CellField(x->x,Γ)⋅nΓ)dΓ))
end

# Every (D-2)-face of the glued interface is shared by two pieces
function is_closed(Γ)
  topo = get_grid_topology(get_active_model(Γ))
  Df = num_cell_dims(Γ)
  all(f -> length(f) == 2, get_faces(topo,Df-1,Df))
end

≈ₘ(a,b) = isapprox(a,b;rtol=1e-12,atol=1e-12)

for case in cases()
  name, model, geo, closed = case
  @testset "$name" begin
    ls, rp = LevelSetCutter(), RecursivePlaneCutter(atol=get(case,5,1e-12))

    @test compute_bgcell_to_inoutcut(rp,model,geo) == compute_bgcell_to_inoutcut(ls,model,geo)
    @test compute_bgfacet_to_inoutcut(rp,model,geo) == compute_bgfacet_to_inoutcut(ls,model,geo)

    cut_ls, cut_rp = cut(ls,model,geo), cut(rp,model,geo)
    @test cut_rp.oid_to_ls == cut_ls.oid_to_ls
    @test cut_rp.ls_to_bgcell_to_inoutcut == cut_ls.ls_to_bgcell_to_inoutcut

    for io in (PHYSICAL_IN,PHYSICAL_OUT)
      @test sum(∫(1)Measure(Triangulation(cut_rp,io),2)) ≈ₘ sum(∫(1)Measure(Triangulation(cut_ls,io),2))
    end
    for io in (ACTIVE_IN,ACTIVE_OUT)
      @test num_cells(Triangulation(cut_rp,io)) == num_cells(Triangulation(cut_ls,io))
    end
    @test num_cells(GhostSkeleton(cut_rp)) == num_cells(GhostSkeleton(cut_ls))

    # Γ of the whole geometry and of every level set
    q_rp, q_ls = boundary_quantities(EmbeddedBoundary(cut_rp)), boundary_quantities(EmbeddedBoundary(cut_ls))
    @test all(map(≈ₘ,q_rp,q_ls))
    for leaf in leaf_names(geo)
      q_rp = boundary_quantities(EmbeddedBoundary(cut_rp,leaf))
      q_ls = boundary_quantities(EmbeddedBoundary(cut_ls,leaf))
      @test all(map(≈ₘ,q_rp,q_ls))
    end
    closed && @test is_closed(EmbeddedBoundary(cut_rp))

    # Cut background facets
    cutf_ls, cutf_rp = cut_facets(ls,model,geo), cut_facets(rp,model,geo)
    @test cutf_rp.ls_to_facet_to_inoutcut == cutf_ls.ls_to_facet_to_inoutcut
    for io in (PHYSICAL_IN,PHYSICAL_OUT)
      @test sum(∫(1)Measure(SkeletonTriangulation(cutf_rp,io),2)) ≈ₘ sum(∫(1)Measure(SkeletonTriangulation(cutf_ls,io),2))
      @test sum(∫(1)Measure(BoundaryTriangulation(cutf_rp,io),2)) ≈ₘ sum(∫(1)Measure(BoundaryTriangulation(cutf_ls,io),2))
    end
  end
end

# Values within atol of zero are zeros (IN), for the cut and the background labels alike
model = CartesianDiscreteModel((0,1,0,1),(4,4))
geo = AnalyticalGeometry(x -> x[1]-0.5+1e-14,name="line")
rp = RecursivePlaneCutter()
cut_rp = cut(rp,model,geo)
@test only(cut_rp.ls_to_bgcell_to_inoutcut) == compute_bgcell_to_inoutcut(rp,model,geo)
@test sum(∫(1)Measure(Triangulation(cut_rp,PHYSICAL_IN),1)) ≈ 0.5
@test sum(∫(1)Measure(EmbeddedBoundary(cut_rp),1)) ≈ 1
labels_ls = compute_bgcell_to_inoutcut(LevelSetCutter(),model,geo)
@test compute_bgcell_to_inoutcut(RecursivePlaneCutter(atol=0),model,geo) == labels_ls
@test compute_bgcell_to_inoutcut(rp,model,geo) != labels_ls # the nodes on x = 0.5 are IN

end # module
