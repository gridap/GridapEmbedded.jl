
using Gridap
using GridapEmbedded
using Test
using LinearAlgebra

include("GridapEmbeddedFixes.jl")

function coarse_active_triangulation(mH,mh,mh_cut)
  # Build the coarse active triangulation 
  # A coarse cell is active if it has at least one active child
  mask_coarse_cells = Vector{Bool}(undef,num_cells(mH))
  mask_coarse_cells .= false
  n2o = mh.glue.n2o_faces_map
  inoutcut = mh_cut.ls_to_bgcell_to_inoutcut[1]
  for fine_cell in 1:length(inoutcut)
    if inoutcut[fine_cell] in (GridapEmbedded.Interfaces.CUT,GridapEmbedded.Interfaces.IN)
      mask_coarse_cells[n2o[end][fine_cell]] = true
    end
  end
  Triangulation(mH,mask_coarse_cells)
end



const R = 0.4

# Weak form Nitsche + ghost penalty (CutFEM paper Sect. 6.1)
const γd = 10.0

function solve_with_agfem(n, order, scheme)
  @assert scheme == :galerkin || scheme == :petrov_galerkin

  # Manufactured solution
  u(x) = x[1] + x[2]^order
  ∇u(x) = ∇(u)(x)
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  # Select geometry
  geom = disk(R, x0=Point(0.5,0.5))
  partition = (n,n)

  # Setup background model
  box = get_metadata(geom)
  mH = CartesianDiscreteModel((0,1,0,1),partition)
  dp = box.pmax - box.pmin
  h = dp[1]/n

  mh = Gridap.Adaptivity.refine(mH,order)

  # Cut the background model with the fine mesh
  cutdisc = cut(mh,geom)

  # Cut the background model with the coarse mesh and build aggregates
  coarse_cutdisk = cut(mH,geom)
  strategy = AggregateAllCutCells()
  coarse_aggregates = aggregate(strategy,coarse_cutdisk)

  # Setup integration meshes
  Ωhact = Triangulation(cutdisc,ACTIVE)
  ΩHact = coarse_active_triangulation(mH,mh,cutdisc)
  Ωh = Triangulation(cutdisc,PHYSICAL)
  Γh = EmbeddedBoundary(cutdisc)

  @test ΩHact.tface_to_mface == findall(x->x in (GridapEmbedded.Interfaces.CUT,GridapEmbedded.Interfaces.IN),
                                  coarse_cutdisk.ls_to_bgcell_to_inoutcut[1])

  # Setup normal vectors
  n_Γ = get_normal_vector(Γh)

  # Setup Lebesgue measures
  degree = 2*order+1
  dΩ = Measure(Ωh,degree)
  dΓ = Measure(Γh,degree)

  # Setup FESpace
  trialRefFE=ReferenceFE(lagrangian,Float64,order)
  Vstd = TestFESpace(ΩHact,trialRefFE,conformity=:H1)
  Ustd = TrialFESpace(Vstd)
  Vagg = AgFEMSpace(Ustd,coarse_aggregates)
  Uagg = TrialFESpace(Vagg)

  if (scheme == :petrov_galerkin)
    testRefFE=Gridap.HQkIsoQ1(Float64,num_cell_dims(mH),order)
    Vstd = TestFESpace(ΩHact,testRefFE,conformity=:H1)
    Vagg=Gridap.FESpaces.FESpaceWithLinearConstraints(
      Vstd,
      Uagg.n_fdofs,
      Uagg.n_fmdofs,
      Uagg.mDOF_to_DOF,
      Uagg.DOF_to_mDOFs,
      Uagg.DOF_to_coeffs,
      Uagg.cell_to_lmdof_to_mdof,
      Uagg.cell_to_ldof_to_dof)
  end

  function a(u,v)
    gvh = Gridap.CellData.change_domain(∇(v),ReferenceDomain(), Ωhact, ReferenceDomain())  
    guh = Gridap.CellData.change_domain(∇(u),ReferenceDomain(), Ωhact, ReferenceDomain())
    uh = Gridap.CellData.change_domain(u,ReferenceDomain(), Ωhact, ReferenceDomain())
    vh = Gridap.CellData.change_domain(v,ReferenceDomain(), Ωhact, ReferenceDomain())
    ∫( gvh⋅guh ) * dΩ +
    ∫( (γd/h)*vh*uh  - vh*(n_Γ⋅guh) - (n_Γ⋅gvh)*uh ) * dΓ
  end 

  function l(v)
    gvh = Gridap.CellData.change_domain(∇(v),ReferenceDomain(), Ωhact, ReferenceDomain()) 
    vh = Gridap.CellData.change_domain(v,ReferenceDomain(), Ωhact, ReferenceDomain())
    ∫( vh*f ) * dΩ +
      ∫( (γd/h)*vh*ud - (n_Γ⋅gvh)*ud ) * dΓ
  end

  # FE problem
  op = AffineFEOperator(a,l,Uagg,Vagg)
  # println(cond(Array(op.op.matrix)))
  uh = Gridap.CellData.change_domain(solve(op),ReferenceDomain(),Ωhact,ReferenceDomain())

  e = Gridap.CellData.change_domain(u - uh,ReferenceDomain(),Ωhact,ReferenceDomain())

  # Postprocess
  l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
  h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

  el2 = l2(e)
  eh1 = h1(e)
  ul2 = l2(uh)
  uh1 = h1(uh)

  # writevtk(Ω,"results",cellfields=["uh"=>uh])
  @test el2/ul2 < 1.e-6
  @test eh1/uh1 < 1.e-6
end 

solve_with_agfem(3,1,:galerkin)
solve_with_agfem(3,2,:galerkin)
solve_with_agfem(3,4,:galerkin)
solve_with_agfem(3,8,:galerkin) # Fails

solve_with_agfem(3,1,:petrov_galerkin)
solve_with_agfem(3,2,:petrov_galerkin)
solve_with_agfem(3,4,:petrov_galerkin)
solve_with_agfem(3,8,:petrov_galerkin) # Fails



# DEBUG statements
# dv = get_fe_basis(V)
# du = get_trial_fe_basis(U)
# gdv = ∇(dv)
# gdu = ∇(du)
# gvh = Gridap.CellData.change_domain(gdv,ReferenceDomain(), Ωhact, ReferenceDomain())  
# guh = Gridap.CellData.change_domain(gdu,ReferenceDomain(), Ωhact, ReferenceDomain())
# gvh_guh=gvh⋅guh
# x=Gridap.CellData.get_cell_points(dΩ.quad)
# gvh_Ω = Gridap.CellData.change_domain(gvh,ReferenceDomain(), Ωh, ReferenceDomain())
# guh_Ω = Gridap.CellData.change_domain(gvh,ReferenceDomain(), Ωh, ReferenceDomain())
# gx=gvh(x)
# gy=guh(x)
# gvh_guh(x)
# gx.b
# gy.b.args[1]
# vh = Gridap.CellData.change_domain(dv,ReferenceDomain(), Ωhact, ReferenceDomain())
# fh = CellField(f,Ωhact)
# (f*vh)(x)
# f_vh=f*vh
# f_vh_Ωh=Gridap.CellData.change_domain(f_vh.args[1],ReferenceDomain(), Ωh, ReferenceDomain()) 
# gx = f_vh.args[1](x)
# fx = f_vh.args[2](x)
# f_vh(x)
# a(du,dv)
# l(dv) 