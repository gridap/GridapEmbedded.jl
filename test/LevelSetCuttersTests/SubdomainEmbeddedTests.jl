module SubdomainEmbeddedTests

using Test
using LinearAlgebra
using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.MultiField
using GridapEmbedded
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.LevelSetCutters: DifferentiableTriangulation

function nudge_zeros!(φh)
  x_φ = get_free_dof_values(φh)
  idx = findall(isapprox(0.0; atol=1e-10), x_φ)
  x_φ[idx] .+= 100 * eps(eltype(x_φ))
  return φh
end

function check_directional_derivative(f, x, g_ad; ε=1e-6, tol=1e-6, ndirs=3)
  for k in 1:ndirs
    d = [sin(k * 1.7 * i + 0.3 * k) for i in eachindex(x)]
    d ./= norm(d)
    fd = (f(x .+ ε .* d) - f(x .- ε .* d)) / (2 * ε)
    @test isapprox(dot(g_ad, d), fd; atol=tol, rtol=tol)
  end
end

############################################################################################
# 1. Subdomain cutting basics, CSG, GhostSkeleton, merge_triangulations, dimension checks
#    (tested for both LevelSetCutter and RecursivePlaneCutter)
############################################################################################
@testset "Subdomain basics & CSG (2D) - $cutter" for cutter in (LevelSetCutter(), RecursivePlaneCutter())
  model = simplexify(CartesianDiscreteModel((0, 1, 0, 1), (10, 10)))
  sq_ids_quads = [34:37; 44:47; 54:57; 64:67]
  sq_ids = vcat([2*i-1:2*i for i in sq_ids_quads]...)
  remaining_ids = setdiff(1:num_cells(model), sq_ids)
  Ω_φ = Triangulation(model, sq_ids)
  Ω_r = Triangulation(model, remaining_ids)

  V_φ = TestFESpace(Ω_φ, ReferenceFE(lagrangian, Float64, 1))
  φh1 = nudge_zeros!(interpolate(x -> x[1] - x[2] - 0.13, V_φ))
  φh2 = nudge_zeros!(interpolate(x -> sqrt((x[1] - 0.5)^2 + (x[2] - 0.5)^2) - 0.17, V_φ))

  geo1 = DiscreteGeometry(φh1, Ω_φ)
  geo2 = DiscreteGeometry(φh2, Ω_φ)
  geo_csg = union(geo1, geo2)
  geo_diff = setdiff(geo1, geo2)
  geo_comp = !geo1
  geo_ana = disk(0.17, x0=Point(0.5, 0.5))

  # Dimension check: cutting a lower-dimensional triangulation must fail
  Γ_bdry = BoundaryTriangulation(model)
  @test_throws AssertionError cut(cutter, Γ_bdry, geo_ana)
  @test_throws AssertionError cut_facets(cutter, Γ_bdry, geo_ana)

  for geo in (geo1, geo_csg, geo_diff, geo_comp, geo_ana)
    cutgeo = cut(cutter, Ω_φ, geo)
    cutfacet = cut_facets(cutter, Ω_φ, geo)

    @test compute_bgcell_to_inoutcut(cutter, Ω_φ, geo) == compute_bgcell_to_inoutcut(cutgeo, geo)
    @test compute_bgfacet_to_inoutcut(cutter, Ω_φ, geo) == compute_bgfacet_to_inoutcut(cutfacet, geo)

    Ω_φ_act = Triangulation(cutgeo, ACTIVE)
    Ω_φ_phys = Triangulation(cutgeo, PHYSICAL)
    Ω_φ_out = Triangulation(cutgeo, PHYSICAL_OUT)
    Γ = EmbeddedBoundary(cutgeo)
    Γg = GhostSkeleton(cutgeo)

    @test get_background_model(Ω_φ_act) === model
    @test get_background_model(Ω_φ_phys) === model
    @test get_background_model(Γ) === model
    @test get_background_model(Γg) === model

    # Volume conservation strictly within Ω_φ
    vol_in = sum(∫(1) * Measure(Ω_φ_phys, 2))
    vol_out = sum(∫(1) * Measure(Ω_φ_out, 2))
    vol_φ = sum(∫(1) * Measure(Ω_φ, 2))
    @test vol_in + vol_out ≈ vol_φ atol=1e-12

    # Active cells subset of sq_ids; GhostSkeleton strictly interior to Ω_φ
    act_cells = get_glue(Ω_φ_act, Val(2)).tface_to_mface
    @test issubset(act_cells, sq_ids)
    if num_cells(Γg) > 0
      @test issubset(get_glue(Γg.⁺, Val(2)).tface_to_mface, sq_ids)
      @test issubset(get_glue(Γg.⁻, Val(2)).tface_to_mface, sq_ids)
    end

    # merge_triangulations builds an injective BodyFittedTriangulation and rejects overlaps
    Ωact = merge_triangulations(Ω_φ_act, Ω_r)
    @test get_background_model(Ωact) === model
    @test sort(collect(get_glue(Ωact, Val(2)).tface_to_mface)) == sort(union(act_cells, remaining_ids))
    @test !isnothing(get_glue(Ωact, Val(2)).mface_to_tface)
    if num_cells(Ω_φ_act) > 0
      @test_throws ErrorException merge_triangulations(Ω_φ_act, Ω_φ)
    end

    # Cut facets via Interface(Ω_r, Ω_φ) and Boundary(Ω_r)
    Λ = Interface(Ω_r, Ω_φ)
    Γ_int_in = Boundary(Λ.⁺, cutfacet, PHYSICAL_IN)
    Γ_int_out = Boundary(Λ.⁺, cutfacet, PHYSICAL_OUT)
    len_in = num_cells(Γ_int_in) > 0 ? sum(∫(1) * Measure(Γ_int_in, 2)) : 0.0
    len_out = num_cells(Γ_int_out) > 0 ? sum(∫(1) * Measure(Γ_int_out, 2)) : 0.0
    len_tot = sum(∫(1) * Measure(Λ.⁺, 2))
    @test len_in + len_out ≈ len_tot atol=1e-12

    Γ_bdry_r_out = Boundary(Boundary(Ω_r), cutfacet, PHYSICAL_OUT)
    len_out_r = num_cells(Γ_bdry_r_out) > 0 ? sum(∫(1) * Measure(Γ_bdry_r_out, 2)) : 0.0
    @test len_out_r ≈ len_out atol=1e-12
  end

  # Cross-check RecursivePlaneCutter vs LevelSetCutter and full-model DiscreteGeometry on subdomain
  V_full = TestFESpace(model, ReferenceFE(lagrangian, Float64, 1))
  φh_full = nudge_zeros!(interpolate(x -> x[1] - x[2] - 0.13, V_full))
  geo_full = DiscreteGeometry(φh_full, model)
  cutgeo_full = cut(cutter, Ω_φ, geo_full)
  cutfacet_full = cut_facets(cutter, cutgeo_full, geo_full)
  @test sum(∫(1) * Measure(Triangulation(cutgeo_full, PHYSICAL), 2)) ≈
        sum(∫(1) * Measure(Triangulation(cut(cutter, Ω_φ, geo1), PHYSICAL), 2)) atol=1e-12
  @test sum(∫(1) * Measure(Boundary(Interface(Ω_r, Ω_φ).⁺, cutfacet_full, PHYSICAL_OUT), 2)) ≈
        sum(∫(1) * Measure(Boundary(Interface(Ω_r, Ω_φ).⁺, cut_facets(cutter, Ω_φ, geo1), PHYSICAL_OUT), 2)) atol=1e-12
end

############################################################################################
# 2. Subdomain CutFEM Poisson PDE (2D) for both cutters + Full AD verification (LevelSetCutter)
############################################################################################
@testset "Subdomain CutFEM Poisson (2D) - $cutter" for cutter in (LevelSetCutter(), RecursivePlaneCutter())
  n = 10
  h = 1.0 / n
  model = simplexify(CartesianDiscreteModel((0, 1, 0, 1), (n, n)))

  sq_ids_quads = [34:37; 44:47; 54:57; 64:67]
  sq_ids = vcat([2*i-1:2*i for i in sq_ids_quads]...)
  remaining_ids = setdiff(1:num_cells(model), sq_ids)
  Ω_φ = Triangulation(model, sq_ids)
  Ω_r = Triangulation(model, remaining_ids)

  V_φ = TestFESpace(Ω_φ, ReferenceFE(lagrangian, Float64, 1))
  φh = nudge_zeros!(interpolate(((x, y),) -> x - y - 0.13, V_φ))

  geo = DiscreteGeometry(φh, Ω_φ)
  cutgeo = cut(cutter, Ω_φ, geo)
  cutfacet = cut_facets(cutter, Ω_φ, geo)

  u(x) = x[1] + 2 * x[2] - 0.5
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  Ω_φ_act = Triangulation(cutgeo, ACTIVE)
  Ω_φ_phys = Triangulation(cutgeo, PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)
  Λ_rφ = Interface(Ω_r, Ω_φ)
  Γ_∂Ωr = Boundary(Λ_rφ.⁺, cutfacet, PHYSICAL_OUT)
  Γg = GhostSkeleton(cutgeo)

  n_Γ = get_normal_vector(Γ)
  n_Γg = get_normal_vector(Γg)
  n_∂Ωr = get_normal_vector(Γ_∂Ωr)

  order = 1
  degree = 2 * order
  Ω_phys = merge_triangulations(Ω_φ_phys, Ω_r)
  Ω_phys_rev = merge_triangulations(Ω_r, Ω_φ_phys)
  @test Ω_phys isa Gridap.Geometry.AppendedTriangulation
  @test Ω_phys.a isa SubFaceTriangulation
  @test Ω_phys.b isa Gridap.Geometry.BodyFittedTriangulation
  @test Ω_phys_rev.a isa SubFaceTriangulation
  @test Ω_phys_rev.b isa Gridap.Geometry.BodyFittedTriangulation

  dΩ_phys = Measure(Ω_phys, degree)
  dΩ_φ_phys = Measure(Ω_φ_phys, degree)
  dΩ_r = Measure(Ω_r, degree)
  @test sum(∫(1) * dΩ_phys) ≈ sum(∫(1) * dΩ_φ_phys + ∫(1) * dΩ_r) atol=1e-12
  @test sum(∫(1) * Measure(Ω_phys_rev, degree)) ≈ sum(∫(1) * dΩ_phys) atol=1e-12
  dΓg = Measure(Γg, degree)
  quad = Quadrature(duffy, degree)
  dΓ = Measure(Γ, quad)
  dΓ_∂Ωr = Measure(Γ_∂Ωr, degree)

  Ωact = merge_triangulations(Ω_φ_act, Ω_r)
  V = TestFESpace(Ωact, ReferenceFE(lagrangian, Float64, order), conformity=:H1, dirichlet_tags="boundary")
  U = TrialFESpace(V, ud)

  γd = 10.0
  γg = 0.1

  a(u, v) =
    ∫(∇(v) ⋅ ∇(u)) * dΩ_phys +
    ∫((γd / h) * v * u - v * (n_Γ ⋅ ∇(u)) - (n_Γ ⋅ ∇(v)) * u) * dΓ +
    ∫((γd / h) * v * u - v * (n_∂Ωr ⋅ ∇(u)) - (n_∂Ωr ⋅ ∇(v)) * u) * dΓ_∂Ωr +
    ∫((γg * h) * jump(n_Γg ⋅ ∇(v)) * jump(n_Γg ⋅ ∇(u))) * dΓg

  l(v) =
    ∫(v * f) * dΩ_phys +
    ∫((γd / h) * v * ud - (n_Γ ⋅ ∇(v)) * ud) * dΓ +
    ∫((γd / h) * v * ud - (n_∂Ωr ⋅ ∇(v)) * ud) * dΓ_∂Ωr

  uh = solve(AffineFEOperator(a, l, U, V))

  e = u - uh
  el2 = sqrt(sum(∫(e * e) * dΩ_phys))
  eh1 = sqrt(sum(∫(e * e + ∇(e) ⋅ ∇(e)) * dΩ_phys))
  @test el2 < 1e-10
  @test eh1 < 1e-9

  # AD verification (LevelSetCutter only)
  if cutter isa LevelSetCutter
    wh = interpolate(x -> sin(π * x[1]) * cos(π * x[2]) + x[1]^2 + 0.5 * x[2]^2, V)

    diff_Ω_phys = DifferentiableTriangulation(Ω_phys, V_φ)
    diff_Γ = DifferentiableTriangulation(Γ, V_φ)
    diff_Γ_∂Ωr = DifferentiableTriangulation(Γ_∂Ωr, V_φ)

    ddiff_Ω_phys = Measure(diff_Ω_phys, degree)
    ddiff_Γ = Measure(diff_Γ, quad)
    ddiff_Γ_∂Ωr = Measure(diff_Γ_∂Ωr, degree)

    # (a) Bulk term on diff_Ω_phys
    J_bulk(φ) = ∫(∇(wh) ⋅ ∇(wh) + wh * wh) * ddiff_Ω_phys
    g_bulk_ad = assemble_vector(gradient(J_bulk, φh), V_φ)
    check_directional_derivative(get_free_dof_values(φh), g_bulk_ad; tol=1e-6) do φ_vec
      cg = cut(Ω_φ, DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ))
      Ωp_i = merge_triangulations(Triangulation(cg, PHYSICAL), Ω_r)
      sum(∫(∇(wh) ⋅ ∇(wh) + wh * wh) * Measure(Ωp_i, degree))
    end

    # (b) Embedded boundary term on diff_Γ (including moving normal)
    function J_gamma(φ)
      n = get_normal_vector(diff_Γ)
      ∫((γd / h) * wh * wh - 2 * wh * (n ⋅ ∇(wh))) * ddiff_Γ
    end
    g_gamma_ad = assemble_vector(gradient(J_gamma, φh), V_φ)
    check_directional_derivative(get_free_dof_values(φh), g_gamma_ad; tol=1e-5) do φ_vec
      cg = cut(Ω_φ, DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ))
      Γ_i = EmbeddedBoundary(cg)
      n_i = get_normal_vector(Γ_i)
      sum(∫((γd / h) * wh * wh - 2 * wh * (n_i ⋅ ∇(wh))) * Measure(Γ_i, quad))
    end

    # (c) Cut background facet term on diff_Γ_∂Ωr (SubFacetBoundaryTriangulation glued to Ω_r)
    function J_facet(φ)
      n = get_normal_vector(diff_Γ_∂Ωr)
      ∫((γd / h) * wh * wh - 2 * wh * (n ⋅ ∇(wh))) * ddiff_Γ_∂Ωr
    end
    g_facet_ad = assemble_vector(gradient(J_facet, φh), V_φ)
    check_directional_derivative(get_free_dof_values(φh), g_facet_ad; tol=1e-5) do φ_vec
      cf = cut_facets(Ω_φ, DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ))
      Γf_i = Boundary(Λ_rφ.⁺, cf, PHYSICAL_OUT)
      nf_i = get_normal_vector(Γf_i)
      sum(∫((γd / h) * wh * wh - 2 * wh * (nf_i ⋅ ∇(wh))) * Measure(Γf_i, degree))
    end

    # (d) Full combined functional j(u, v, φ)
    function j_full(u_f, v_f, φ)
      n_Γ_ad = get_normal_vector(diff_Γ)
      n_∂Ωr_ad = get_normal_vector(diff_Γ_∂Ωr)
      ∫(∇(v_f) ⋅ ∇(u_f)) * ddiff_Ω_phys +
      ∫((γd / h) * v_f * u_f - v_f * (n_Γ_ad ⋅ ∇(u_f)) - (n_Γ_ad ⋅ ∇(v_f)) * u_f) * ddiff_Γ +
      ∫((γd / h) * v_f * u_f - v_f * (n_∂Ωr_ad ⋅ ∇(u_f)) - (n_∂Ωr_ad ⋅ ∇(v_f)) * u_f) * ddiff_Γ_∂Ωr +
      ∫((γg * h) * jump(n_Γg ⋅ ∇(v_f)) * jump(n_Γg ⋅ ∇(u_f))) * dΓg
    end
    g_full_ad = assemble_vector(gradient(φ -> j_full(wh, wh, φ), φh), V_φ)
    check_directional_derivative(get_free_dof_values(φh), g_full_ad; tol=1e-5) do φ_vec
      g_i = DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ)
      cg = cut(Ω_φ, g_i)
      cf = cut_facets(Ω_φ, g_i)
      Ωp_i = merge_triangulations(Triangulation(cg, PHYSICAL), Ω_r)
      Γ_i = EmbeddedBoundary(cg)
      Γf_i = Boundary(Λ_rφ.⁺, cf, PHYSICAL_OUT)
      n_Γ_i = get_normal_vector(Γ_i)
      n_f_i = get_normal_vector(Γf_i)
      sum(
        ∫(∇(wh) ⋅ ∇(wh)) * Measure(Ωp_i, degree) +
        ∫((γd / h) * wh * wh - 2 * wh * (n_Γ_i ⋅ ∇(wh))) * Measure(Γ_i, quad) +
        ∫((γd / h) * wh * wh - 2 * wh * (n_f_i ⋅ ∇(wh))) * Measure(Γf_i, degree) +
        ∫((γg * h) * jump(n_Γg ⋅ ∇(wh)) * jump(n_Γg ⋅ ∇(wh))) * dΓg
      )
    end

    # (e) MultiField AD with V_φ on Ω_φ and V on Ωact
    Y = MultiFieldFESpace([V_φ, V])
    φwh = interpolate([φh, wh], Y)
    g_mf = assemble_vector(gradient(((φ, w),) -> j_full(w, w, φ), φwh), Y)
    g_w_ad = assemble_vector(gradient(w -> j_full(w, w, φh), wh), V)
    n_φ = num_free_dofs(V_φ)
    @test norm(g_mf[1:n_φ] - g_full_ad, Inf) < 1e-12
    @test norm(g_mf[n_φ+1:end] - g_w_ad, Inf) < 1e-12
  end
end

############################################################################################
# 3. 3D Subdomain CutFEM (both cutters) + SubFacetBoundaryTriangulation AD (LevelSetCutter)
############################################################################################
@testset "Subdomain CutFEM (3D) - $cutter" for cutter in (LevelSetCutter(), RecursivePlaneCutter())
  n = 4
  h = 1.0 / n
  model = simplexify(CartesianDiscreteModel((0, 1, 0, 1, 0, 1), (n, n, n)))

  cell_coords = get_cell_coordinates(model)
  sq_ids = findall(cell_coords) do pts
    c = sum(pts) / length(pts)
    all(0.25 .<= Tuple(c) .<= 0.75)
  end
  remaining_ids = setdiff(1:num_cells(model), sq_ids)
  Ω_φ = Triangulation(model, sq_ids)
  Ω_r = Triangulation(model, remaining_ids)

  V_φ = TestFESpace(Ω_φ, ReferenceFE(lagrangian, Float64, 1))
  φh = nudge_zeros!(interpolate(x -> x[1] - 0.5 * x[2] + 0.3 * x[3] - 0.41, V_φ))

  geo = DiscreteGeometry(φh, Ω_φ)
  cutgeo = cut(cutter, Ω_φ, geo)
  cutfacet = cut_facets(cutter, Ω_φ, geo)

  Ω_φ_act = Triangulation(cutgeo, ACTIVE)
  Ω_φ_phys = Triangulation(cutgeo, PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)
  Λ_rφ = Interface(Ω_r, Ω_φ)
  Γ_∂Ωr = Boundary(Λ_rφ.⁺, cutfacet, PHYSICAL_OUT)
  Γg = GhostSkeleton(cutgeo)

  Ωact = merge_triangulations(Ω_φ_act, Ω_r)
  Ω_phys = merge_triangulations(Ω_φ_phys, Ω_r)
  order = 1
  degree = 2 * order
  u(x) = x[1] + 2 * x[2] - 3 * x[3] + 0.25
  f(x) = -Δ(u)(x)
  V = TestFESpace(Ωact, ReferenceFE(lagrangian, Float64, order), conformity=:H1, dirichlet_tags="boundary")
  U = TrialFESpace(V, u)

  n_Γ = get_normal_vector(Γ)
  n_Γg = get_normal_vector(Γg)
  n_∂Ωr = get_normal_vector(Γ_∂Ωr)

  dΩ_phys = Measure(Ω_phys, degree)
  dΓ = Measure(Γ, Quadrature(duffy, degree))
  dΓ_∂Ωr = Measure(Γ_∂Ωr, degree)
  dΓg = Measure(Γg, degree)

  γd = 10.0
  γg = 0.1
  a(u, v) =
    ∫(∇(v) ⋅ ∇(u)) * dΩ_phys +
    ∫((γd / h) * v * u - v * (n_Γ ⋅ ∇(u)) - (n_Γ ⋅ ∇(v)) * u) * dΓ +
    ∫((γd / h) * v * u - v * (n_∂Ωr ⋅ ∇(u)) - (n_∂Ωr ⋅ ∇(v)) * u) * dΓ_∂Ωr +
    ∫((γg * h) * jump(n_Γg ⋅ ∇(v)) * jump(n_Γg ⋅ ∇(u))) * dΓg

  l(v) =
    ∫(v * f) * dΩ_phys +
    ∫((γd / h) * v * u - (n_Γ ⋅ ∇(v)) * u) * dΓ +
    ∫((γd / h) * v * u - (n_∂Ωr ⋅ ∇(v)) * u) * dΓ_∂Ωr

  uh = solve(AffineFEOperator(a, l, U, V))
  e = u - uh
  el2 = sqrt(sum(∫(e * e) * dΩ_phys))
  eh1 = sqrt(sum(∫(e * e + ∇(e) ⋅ ∇(e)) * dΩ_phys))
  @test el2 < 1e-10
  @test eh1 < 1e-9

  if cutter isa LevelSetCutter
    wh = interpolate(x -> x[1]^2 + 2 * x[2]^2 - x[3]^2 + x[1] * x[2] * x[3], V)
    diff_Ω_phys = DifferentiableTriangulation(Ω_phys, V_φ)
    diff_Γ = DifferentiableTriangulation(Γ, V_φ)
    diff_Γ_∂Ωr = DifferentiableTriangulation(Γ_∂Ωr, V_φ)

    ddiff_Ω_phys = Measure(diff_Ω_phys, degree)
    ddiff_Γ = Measure(diff_Γ, Quadrature(duffy, degree))
    ddiff_Γ_∂Ωr = Measure(diff_Γ_∂Ωr, degree)

    function J_3d(φ)
      n1 = get_normal_vector(diff_Γ)
      n2 = get_normal_vector(diff_Γ_∂Ωr)
      ∫(wh * wh + ∇(wh) ⋅ ∇(wh)) * ddiff_Ω_phys +
      ∫(wh * wh - wh * (n1 ⋅ ∇(wh))) * ddiff_Γ +
      ∫(wh * wh - wh * (n2 ⋅ ∇(wh))) * ddiff_Γ_∂Ωr
    end
    g_3d_ad = assemble_vector(gradient(J_3d, φh), V_φ)
    check_directional_derivative(get_free_dof_values(φh), g_3d_ad; tol=1e-5) do φ_vec
      g_i = DiscreteGeometry(FEFunction(V_φ, φ_vec), Ω_φ)
      cg = cut(Ω_φ, g_i)
      cf = cut_facets(Ω_φ, g_i)
      Ωp_i = merge_triangulations(Triangulation(cg, PHYSICAL), Ω_r)
      Γ_i = EmbeddedBoundary(cg)
      Γf_i = Boundary(Λ_rφ.⁺, cf, PHYSICAL_OUT)
      n1 = get_normal_vector(Γ_i)
      n2 = get_normal_vector(Γf_i)
      sum(
        ∫(wh * wh + ∇(wh) ⋅ ∇(wh)) * Measure(Ωp_i, degree) +
        ∫(wh * wh - wh * (n1 ⋅ ∇(wh))) * Measure(Γ_i, Quadrature(duffy, degree)) +
        ∫(wh * wh - wh * (n2 ⋅ ∇(wh))) * Measure(Γf_i, degree)
      )
    end
  end
end

end # module
