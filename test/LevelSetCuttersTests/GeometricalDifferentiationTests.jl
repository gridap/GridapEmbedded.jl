module GeometricalDifferentiationTests
############################################################################################
# These tests are meant to verify the correctness differentiation of functionals w.r.t the 
# level set defining the cut domain. 
# They are based on the following work:
# "Level-set topology optimisation with unfitted finite elements and automatic shape differentiation"
#   by Z. J. Wegert, J. Manyer, C. Mallon, S. Badia, V. J. Challis, CMAME (2025)
############################################################################################
using Test, FiniteDiff

using Gridap, Gridap.Geometry, Gridap.Adaptivity, Gridap.Arrays
using GridapEmbedded, GridapEmbedded.LevelSetCutters, GridapEmbedded.Interfaces

using GridapEmbedded.Interfaces: get_conormal_vector
using GridapEmbedded.Interfaces: get_subfacet_normal_vector
using GridapEmbedded.Interfaces: get_ghost_normal_vector

using GridapEmbedded.LevelSetCutters: DifferentiableTriangulation


using Gridap.Fields, Gridap.Polynomials
function Arrays.return_cache(
  fg::Fields.FieldGradientArray{1,Polynomials.MonomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {D,V}
  xi = testitem(x)
  T = gradient_type(V,xi)
  Polynomials._return_cache(fg,x,T,Val(false))
end

function Arrays.evaluate!(
  cache,
  fg::Fields.FieldGradientArray{1,Polynomials.MonomialBasis{D,V}},
  x::AbstractVector{<:Point}) where {D,V}
  Polynomials._evaluate!(cache,fg,x,Val(false))
end

# We general a simplicial model where the simplices are created in a symmetric way using 
# varycentric refinement of QUADs and HEXs.
function generate_model(D,n)
  domain = (D==2) ? (0,1,0,1) : (0,1,0,1,0,1)
  cell_partition = (D==2) ? (n,n) : (n,n,n)
  base_model = UnstructuredDiscreteModel((CartesianDiscreteModel(domain,cell_partition)))
  ref_model = refine(base_model, refinement_method = "barycentric")
  model = ref_model.model
  return model
end

function level_set(shape::Symbol;N=4)
  if shape == :square
    x -> max(abs(x[1]-0.5),abs(x[2]-0.5))-0.25 # Square
  elseif shape == :corner_2d
    x -> ((x[1]-0.5)^N+(x[2]-0.5)^N)^(1/N)-0.25 # Curved corner
  elseif shape == :diamond
    x -> abs(x[1]-0.5)+abs(x[2]-0.5)-0.25-0/n/10 # Diamond
  elseif shape == :circle
    x -> sqrt((x[1]-0.5)^2+(x[2]-0.5)^2)-0.5223 # Circle
  elseif shape == :circle_2
    x -> sqrt((x[1]-0.5)^2+(x[2]-0.5)^2)-0.23 # Circle
  elseif shape == :square_prism
    x -> max(abs(x[1]-0.5),abs(x[2]-0.5),abs(x[3]-0.5))-0.25 # Square prism
  elseif shape == :corner_3d
    x -> ((x[1]-0.5)^N+(x[2]-0.5)^N+(x[3]-0.5)^N)^(1/N)-0.25 # Curved corner
  elseif shape == :diamond_prism
    x -> abs(x[1]-0.5)+abs(x[2]-0.5)+abs(x[3]-0.5)-0.25-0/n/10 # Diamond prism
  elseif shape == :sphere
    x -> sqrt((x[1]-0.5)^2+(x[2]-0.5)^2+(x[3]-0.5)^2)-0.53 # Sphere
  elseif shape == :sphere_2
    x -> sqrt((x[1]-0.5)^2+(x[2]-0.5)^2+(x[3]-0.5)^2)-0.23 # Sphere
  elseif shape == :regular_2d
    x -> cos(2π*x[1])*cos(2π*x[2])-0.11 # "Regular" LSF
  elseif shape == :regular_3d
    x -> cos(2π*x[1])*cos(2π*x[2])*cos(2π*x[3])-0.11 # "Regular" LSF
  else
    error("Unknown shape")
  end
end

function main(
  model,ls::Function,f::Function;
  vtk=false,
  name="embedded",
  verbose=false,
  fdm=false
)
  order = 1
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_φ = TestFESpace(model,reffe)

  U = TestFESpace(model,reffe)

  φh = interpolate(ls,V_φ)
  fh = interpolate(f,V_φ)
  uh = interpolate(x->x[1]+x[2],U)

  # Correction if level set is on top of a node
  x_φ = get_free_dof_values(φh)
  idx = findall(isapprox(0.0;atol=10^-10),x_φ)
  !isempty(idx) && @info "Correcting level values!"
  x_φ[idx] .+= 100*eps(eltype(x_φ))

  geo = DiscreteGeometry(φh,model)
  cutgeo = cut(model,geo)

  # A.1) Volume integral

  Ω = Triangulation(cutgeo,PHYSICAL_IN)
  Ω_AD = DifferentiableTriangulation(Ω,V_φ)
  dΩ = Measure(Ω_AD,2*order)

  Γ = EmbeddedBoundary(cutgeo)
  n_Γ = get_normal_vector(Γ)
  dΓ = Measure(Γ,2*order)

  J_bulk(φ) = ∫(fh)dΩ
  dJ_bulk_AD = gradient(J_bulk,φh)
  dJ_bulk_AD_vec = assemble_vector(dJ_bulk_AD,V_φ)

  dJ_bulk_exact(q) = ∫(-fh*q/(abs(n_Γ ⋅ ∇(φh))))dΓ
  dJ_bulk_exact_vec = assemble_vector(dJ_bulk_exact,V_φ)

  abs_error = norm(dJ_bulk_AD_vec - dJ_bulk_exact_vec,Inf)

  if fdm
    function J_fdm_bulk(φ)
      φh = FEFunction(V_φ,φ)
      cutgeo = cut(model,DiscreteGeometry(φh,model))
      Ω = Triangulation(cutgeo,PHYSICAL_IN)
      dΩ = Measure(Ω,2*order)
      sum(∫(fh)dΩ)
    end
    dJ_FD = FiniteDiff.finite_difference_gradient(J_fdm_bulk,get_free_dof_values(φh))

    abs_error_fdm = norm(dJ_bulk_AD_vec - dJ_FD,Inf)
  end

  if verbose
    println("A.1) Volume integral:")
    println("  - norm(dJ_AD - dJ_exact,Inf) = ",abs_error)
    fdm && println("  - norm(dJ_AD - dJ_FDM,Inf) = ",abs_error_fdm)
  end

  @test abs_error < 1e-10

  # A.1.1) Volume integral with another field

  J_bulk_1(u,φ) = ∫(u+fh)dΩ
  dJ_bulk_1_AD = gradient(φ->J_bulk_1(uh,φ),φh)
  dJ_bulk_1_AD_vec = assemble_vector(dJ_bulk_1_AD,V_φ)

  dJ_bulk_1_exact(q,u) = ∫(-(u+fh)*q/(abs(n_Γ ⋅ ∇(φh))))dΓ
  dJ_bulk_1_exact_vec = assemble_vector(q->dJ_bulk_1_exact(q,uh),V_φ)
  @test norm(dJ_bulk_1_AD_vec - dJ_bulk_1_exact_vec) < 1e-10

  dJ_bulk_1_AD_in_u = gradient(u->J_bulk_1(u,φh),uh)
  dJ_bulk_1_AD_in_u_vec = assemble_vector(dJ_bulk_1_AD_in_u,U)
  dJ_bulk_1_exact_in_u(q,u) = ∫(q)dΩ
  dJ_bulk_1_exact_in_u_vec = assemble_vector(q->dJ_bulk_1_exact_in_u(q,uh),U)

  abs_error = norm(dJ_bulk_1_AD_in_u_vec - dJ_bulk_1_exact_in_u_vec,Inf)
  if verbose
    println("A.1.1) Volume integral with another field:")
    println("  - norm(dJ_AD - dJ_exact,Inf) = ",abs_error)
  end
  @test abs_error < 1e-10

  # A.2) Volume integral

  g(fh) = ∇(fh)⋅∇(fh)
  J_bulk2(φ) = ∫(g(fh))dΩ
  dJ_bulk_AD2 = gradient(J_bulk2,φh)
  dJ_bulk_AD_vec2 = assemble_vector(dJ_bulk_AD2,V_φ)

  dJ_bulk_exact2(q) = ∫(-g(fh)*q/(abs(n_Γ ⋅ ∇(φh))))dΓ
  dJ_bulk_exact_vec2 = assemble_vector(dJ_bulk_exact2,V_φ)

  abs_error = norm(dJ_bulk_AD_vec2 - dJ_bulk_exact_vec2,Inf)

  if verbose
    println("A.2) Volume integral with grad of fields:")
    println("  - norm(dJ_AD - dJ_exact,Inf) = ",abs_error)
  end

  @test abs_error < 1e-10

  # B.1) Facet integral

  Γ = EmbeddedBoundary(cutgeo)
  Γ_AD = DifferentiableTriangulation(Γ,V_φ)
  Λ = Skeleton(Γ)
  Σ = Boundary(Γ)

  dΓ = Measure(Γ,2*order)
  dΛ = Measure(Λ,2*order)
  dΣ = Measure(Σ,2*order)

  n_Γ = get_normal_vector(Γ)

  n_S_Λ = get_normal_vector(Λ)
  m_k_Λ = get_conormal_vector(Λ)
  ∇ˢφ_Λ = Operation(abs)(n_S_Λ ⋅ ∇(φh).plus)

  n_S_Σ = get_normal_vector(Σ)
  m_k_Σ = get_conormal_vector(Σ)
  ∇ˢφ_Σ = Operation(abs)(n_S_Σ ⋅ ∇(φh))

  dΓ_AD = Measure(Γ_AD,2*order)
  J_int(φ) = ∫(fh)dΓ_AD
  dJ_int_AD = gradient(J_int,φh)
  dJ_int_AD_vec = assemble_vector(dJ_int_AD,V_φ)

  dJ_int_exact(w) = ∫((-n_Γ⋅∇(fh))*w/(abs(n_Γ ⋅ ∇(φh))))dΓ +
                    ∫(-n_S_Λ ⋅ (jump(fh*m_k_Λ) * mean(w) / ∇ˢφ_Λ))dΛ +
                    ∫(-n_S_Σ ⋅ (fh*m_k_Σ * w / ∇ˢφ_Σ))dΣ
  dJ_int_exact_vec = assemble_vector(dJ_int_exact,V_φ)

  abs_error = norm(dJ_int_AD_vec - dJ_int_exact_vec,Inf)

  if fdm
    Ω_data = EmbeddedCollection(model,φh) do cutgeo,_,_
      Γ_AD = DifferentiableTriangulation(EmbeddedBoundary(cutgeo),V_φ)
      (;:dΓ_AD => Measure(Γ_AD,2*order))
    end
    function J_fdm_surf(φ)
      sum(∫(fh)Ω_data.dΓ_AD)
    end
    dJ_surf_FD = FiniteDiff.finite_difference_gradient(J_fdm_surf,get_free_dof_values(φh))

    abs_error_fdm = norm(dJ_int_AD_vec - dJ_surf_FD,Inf)
  end

  if verbose
    println("B.1) Surface integral:")
    println("  - norm(dJ_AD - dJ_exact,Inf) = ",abs_error)
    fdm && println("  - norm(dJ_AD - dJ_FDM,Inf) = ",abs_error_fdm)
  end
  @test abs_error < 1e-10

  # B.2) Facet integral
  g(fh) = ∇(fh)⋅∇(fh)
  ∇g(∇∇f,∇f) = ∇∇f⋅∇f + ∇f⋅∇∇f

  J_int2(φ) = ∫(g(fh))dΓ_AD
  dJ_int_AD2 = gradient(J_int2,φh)
  dJ_int_AD_vec2 = assemble_vector(dJ_int_AD2,V_φ)

  dJ_int_exact2(w) = ∫((-n_Γ⋅ (∇g ∘ (∇∇(fh),∇(fh))))*w/(abs(n_Γ ⋅ ∇(φh))))dΓ +
                    ∫(-n_S_Λ ⋅ (jump(g(fh)*m_k_Λ) * mean(w) / ∇ˢφ_Λ))dΛ +
                    ∫(-n_S_Σ ⋅ (g(fh)*m_k_Σ * w / ∇ˢφ_Σ))dΣ
  dJ_int_exact_vec2 = assemble_vector(dJ_int_exact2,V_φ)

  abs_error = norm(dJ_int_AD_vec2 - dJ_int_exact_vec2,Inf)
  if verbose
    println("B.2) Surface integral with grad of other fields:")
    println("  - norm(dJ_AD - dJ_exact,Inf) = ",abs_error)
  end
  @test abs_error < 1e-10

  if vtk
    path = "results/$(name)/"
    mkpath(path)
    Ω_bg = Triangulation(model)
    writevtk(
      Ω_bg,"$(path)results",
      cellfields = [
        "φh" => φh,"∇φh" => ∇(φh),
        "dJ_bulk_AD" => FEFunction(V_φ,dJ_bulk_AD_vec),
        "dJ_bulk_exact" => FEFunction(V_φ,dJ_bulk_exact_vec),
        "dJ_int_AD" => FEFunction(V_φ,dJ_int_AD_vec),
        "dJ_int_exact" => FEFunction(V_φ,dJ_int_exact_vec)
      ],
      celldata = [
        "inoutcut" => GridapEmbedded.Interfaces.compute_bgcell_to_inoutcut(model,geo)
      ];
      append = false
    )

    writevtk(
      Ω, "$(path)omega"; append = false
    )
    writevtk(
      Γ, "$(path)gamma"; append = false
    )

    n_∂Ω_Λ = get_subfacet_normal_vector(Λ)
    n_k_Λ  = get_ghost_normal_vector(Λ)
    writevtk(
      Λ, "$(path)_lambda",
      cellfields = [
        "n_∂Ω.plus" => n_∂Ω_Λ.plus,"n_∂Ω.minus" => n_∂Ω_Λ.minus,
        "n_k.plus" => n_k_Λ.plus,"n_k.minus" => n_k_Λ.minus,
        "n_S" => n_S_Λ,
        "m_k.plus" => m_k_Λ.plus,"m_k.minus" => m_k_Λ.minus,
        "∇ˢφ" => ∇ˢφ_Λ,
        "∇φh_Γs_plus" => ∇(φh).plus,"∇φh_Γs_minus" => ∇(φh).minus,
        "jump(fh*m_k)" => jump(fh*m_k_Λ)
      ];
      append = false
    )

    if num_cells(Σ) > 0
      n_∂Ω_Σ = get_subfacet_normal_vector(Σ)
      n_k_Σ  = get_ghost_normal_vector(Σ)
      writevtk(
        Σ, "$(path)_sigma",
        cellfields = [
          "n_∂Ω" => n_∂Ω_Σ, "n_k" => n_k_Σ,
          "n_S" => n_S_Σ, "m_k" => m_k_Σ,
          "∇ˢφ" => ∇ˢφ_Σ, "∇φh_Γs" => ∇(φh),
        ];
        append = false
      )
    end
  end
end

## Concering integrals of the form `φ->∫(f ⋅ n(φ))dΓ(φ)`
function main_normal(
  model,ls::Function,f_vec::Function;
  vtk=false,
  name="flux integrals",
  run_test=true,
  verbose=false,
  fdm=false
)
  order = 1
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_φ = TestFESpace(model,reffe)

  φh = interpolate(ls,V_φ)

  # Correction if level set is on top of a node
  x_φ = get_free_dof_values(φh)
  idx = findall(isapprox(0.0;atol=10^-10),x_φ)
  !isempty(idx) && @info "Correcting level values!"
  x_φ[idx] .+= 100*eps(eltype(x_φ))

  geo = DiscreteGeometry(φh,model)
  cutgeo = cut(model,geo)

  Γ = EmbeddedBoundary(cutgeo)
  n_Γ = get_normal_vector(Γ)
  Γ_AD = DifferentiableTriangulation(Γ,V_φ)
  dΓ_AD = Measure(Γ_AD,2*order)
  dΓ = Measure(Γ,2*order)

  fh_Γ = CellField(f_vec,Γ)
  fh_Γ_AD = CellField(f_vec,Γ_AD)

  function J_int(φ)
    n = get_normal_vector(Γ_AD)
    ∫(fh_Γ_AD⋅n)dΓ_AD
  end
  dJ_int_AD = gradient(J_int,φh)
  dJ_int_AD_vec = assemble_vector(dJ_int_AD,V_φ)

  _n(∇φ) = ∇φ/(10^-20+norm(∇φ))
  dJ_int_phi = ∇(φ->∫(fh_Γ_AD ⋅ (_n ∘ (∇(φ))))dΓ_AD,φh)
  dJh_int_phi = assemble_vector(dJ_int_phi,V_φ)

  run_test && @test norm(dJ_int_AD_vec - dJh_int_phi) < 1e-10

  # Analytic
  # Note: currently, the analytic result is only valid on closed domains thanks
  #   to the divergence theorem. I think it would take significant work to compute
  #   the analytic derivative generally as we can't rely on divergence theorem to
  #   rewrite it in a convenient way. As a result, we don't have an analytic result
  #   for general cases such as ∫( f(n(φ)) )dΓ(φ), nor the case when Γ intersects
  #   ∂D. Thankfully, we have AD instead ;)
  # Note 2: For the case that Γ does intersect the surface, the result is correct
  #   everywhere except on the intersection.

  fh_Γ = CellField(f_vec,Γ)
  fh_Γ_AD = CellField(f_vec,Γ_AD)

  # Note: this comes from rewriting via the divergence theorem:
  #         ∫(f ⋅ n(φ))dΓ(φ) = ∫(∇⋅f)dΩ(φ)
  dJ_int_exact3(w) = ∫(-(∇⋅(fh_Γ))*w/(abs(n_Γ ⋅ ∇(φh))))dΓ
  dJh_int_exact3 = assemble_vector(dJ_int_exact3,V_φ)

  # Finite diff
  if fdm
    function J_fdm_surf(φ)
      φh = FEFunction(V_φ,φ)
      cutgeo = cut(model,DiscreteGeometry(φh,model))
      Γ = EmbeddedBoundary(cutgeo)
      dΓ = Measure(Γ,2*order)
      n = get_normal_vector(Γ)
      f = CellField(f_vec,Γ)
      n = get_normal_vector(Γ)
      sum(∫(f⋅n)dΓ)
    end
    dJ_FD = FiniteDiff.finite_difference_gradient(J_fdm_surf,get_free_dof_values(φh))

    abs_error_fdm = norm(dJ_int_AD_vec - dJ_FD,Inf)
  end
  abs_error = norm(dJh_int_exact3 - dJ_int_AD_vec,Inf)

  if verbose
    println("C) Flux integral:")
    println("  - norm(dJ_AD - dJ_exact,Inf) = ",abs_error)
    fdm && println("  - norm(dJ_AD - dJ_FDM,Inf) = ",abs_error_fdm)
  end

  run_test && @test abs_error < 1e-10

  if vtk
    path = "results/$(name)/"
    mkpath(path)
    Ω_bg = Triangulation(model)
    writevtk(Ω_bg,path*"Results",cellfields=[
      "dJ_AD"=>FEFunction(V_φ,dJ_int_AD_vec),
      "dJ_AD_with_phi"=>FEFunction(V_φ,dJh_int_phi),
      "dJ_exact"=>FEFunction(V_φ,dJh_int_exact3)
      ])
    writevtk(Γ,path*"Gamma")
  end
end

# Finite difference verification of gradients of integrals of the form `φ->∫(f(n))dΓ(φ)` and Hessian's.
#   Both of these do not currently have rigorous mathematical counterparts so we verify them
#   with finite differences.
function main_fdm_only_verif(model,ls::Function;
    verbose=false,compute_hess=false,fdm=false
)
  order = 1
  reffe = ReferenceFE(lagrangian,Float64,order)
  V_φ = TestFESpace(model,reffe)
  φh = interpolate(ls,V_φ)

  # Correction if level set is on top of a node
  x_φ = get_free_dof_values(φh)
  idx = findall(isapprox(0.0;atol=10^-10),x_φ)
  !isempty(idx) && @info "Correcting level values!"
  x_φ[idx] .+= 100*eps(eltype(x_φ))

  geo = DiscreteGeometry(φh,model)
  cutgeo = cut(model,geo)

  Γ = EmbeddedBoundary(cutgeo)
  Γ_AD = DifferentiableTriangulation(Γ,V_φ)
  dΓ_AD = Measure(Γ_AD,2*order)

  # Sec 6.2 - 10.1007/s00466-017-1383-6
  g((x,y)) = x - 1/10*sin(2π*y/6)
  gh = interpolate(g,V_φ)
  _n_g(∇g) = ∇g/(10^-20+norm(∇g))
  n_g = _n_g ∘ ∇(gh)
  j(x) = norm(x)^2

  function J_int(φ)
    n = get_normal_vector(Γ_AD)
    ∫(j ∘ (n-n_g))dΓ_AD
  end
  dJ_int_AD = gradient(J_int,φh)
  dJ_int_AD_vec = assemble_vector(dJ_int_AD,V_φ)
  hess = hessian(J_int,φh)
  d²J = assemble_matrix(hess,V_φ,V_φ)

  # Finite diff
  if fdm
    function J_fdm_surf(φ)
      φh = FEFunction(V_φ,φ)
      cutgeo = cut(model,DiscreteGeometry(φh,model))
      Γ = EmbeddedBoundary(cutgeo)
      dΓ = Measure(Γ,2*order)
      n = get_normal_vector(Γ)
      sum(∫(j ∘ (n-n_g))dΓ)
    end
    dJ_FD = FiniteDiff.finite_difference_gradient(J_fdm_surf,get_free_dof_values(φh))
    d²J_FD = compute_hess ? FiniteDiff.finite_difference_hessian(J_fdm_surf,get_free_dof_values(φh)) : nothing

    abs_error_fdm = norm(dJ_int_AD_vec - dJ_FD,Inf)
    abs_error_fdm_hess = compute_hess ? norm(d²J - d²J_FD,Inf) : nothing

    if verbose
      println("D) g(n) surf integral:")
      println("  - norm(dJ_AD - dJ_exact,Inf) = ","N/A")
      println("  - norm(dJ_AD - dJ_FDM,Inf) = ",abs_error_fdm)
      compute_hess && println("  - norm(d²J_AD - d²J_FDM,Inf) = ",abs_error_fdm_hess)
    end
  end
end

#######################

# FDM is quite expensive, so we only run if required.

D = 2
n = 10
model = generate_model(D,n)
f(x) = x[1]+x[2]
fvec((x,y)) = VectorValue((1-x)^2,(1-y)^2)
main(model,level_set(:circle_2),f;vtk=false,verbose=true,name="2D_circle/")#,fdm=true)
main(model,level_set(:regular_2d),f;vtk=false,verbose=true,name="2D_reg/")#,fdm=true)
main_normal(model,level_set(:circle_2),fvec;vtk=false,verbose=true,name="2D_circle_flux/",run_test=true)#,fdm=true)
main_normal(model,level_set(:regular_2d),fvec;vtk=false,verbose=true,name="2D_reg_flux/",run_test=false)#,fdm=true) # This will fail as expected
main_fdm_only_verif(model,level_set(:circle_2),verbose=true)#,fdm=true)
main_fdm_only_verif(model,level_set(:regular_2d),verbose=true,compute_hess=true)#,fdm=true)

D = 3
n = 4
model = generate_model(D,n)
φ = level_set(:regular_3d)
f(x) = x[1]+x[2]
fvec2((x,y,z)) = VectorValue((1-x)^2,(1-y)^2,0)
main(model,level_set(:sphere_2),f;vtk=false,verbose=true,name="3D_circle/")#,fdm=true)
main(model,level_set(:regular_3d),f;vtk=false,verbose=true,name="3D_reg/")#,fdm=true)
main_normal(model,level_set(:sphere_2),fvec2;vtk=false,verbose=true,name="3D_circle_flux/",run_test=true)#,fdm=true)
main_normal(model,level_set(:regular_3d),fvec2;vtk=false,verbose=true,name="3D_reg_flux/",run_test=false,)#fdm=true) # This will fail as expected
main_fdm_only_verif(model,level_set(:sphere_2),verbose=true)#,fdm=true)
main_fdm_only_verif(model,level_set(:regular_3d),verbose=true)#,fdm=true)

end