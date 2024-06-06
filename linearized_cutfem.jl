
using Gridap
using GridapEmbedded
using Test

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

# Manufactured solution
u(x) = x[1] + x[2]
∇u(x) = ∇(u)(x)
f(x) = -Δ(u)(x)
ud(x) = u(x)

# Select geometry
const R = 0.2
geom = disk(R, x0=Point(0.5,0.5))
n = 10
partition = (n,n)
order = 1

# Setup background model
box = get_metadata(geom)
mH = CartesianDiscreteModel((0,1,0,1),partition)
dp = box.pmax - box.pmin
const h = dp[1]/n

mh = Gridap.Adaptivity.refine(mH,order)

# Cut the background model with the fine mesh
cutdisc = cut(mh,geom)

# Setup integration meshes
Ωhact = Triangulation(cutdisc,ACTIVE)
ΩHact = coarse_active_triangulation(mH,mh,cutdisc)
Ωh = Triangulation(cutdisc,PHYSICAL)
Γh = EmbeddedBoundary(cutdisc)
Γg = GhostSkeleton(cutdisc)

# Setup normal vectors
n_Γ = get_normal_vector(Γh)
n_Γg = get_normal_vector(Γg)

# Setup Lebesgue measures
degree = 2*order
dΩ = Measure(Ωh,degree)
dΓ = Measure(Γh,degree)
dΓg = Measure(Γg,degree)

# Setup FESpace
trialRefFE=ReferenceFE(lagrangian,Float64,order)
U = TrialFESpace(TestFESpace(ΩHact,trialRefFE,conformity=:H1))
testRefFE=Gridap.HQkIsoQ1(Float64,num_cell_dims(mH),order)
V = TestFESpace(ΩHact,testRefFE,conformity=:H1)

# Weak form Nitsche + ghost penalty (CutFEM paper Sect. 6.1)
const γd = 10.0
const γg = 0.1

function a(u,v)
  gvh = Gridap.CellData.change_domain(∇(v),ReferenceDomain(), Ωhact, ReferenceDomain())  
  guh = Gridap.CellData.change_domain(∇(u),ReferenceDomain(), Ωhact, ReferenceDomain())
  uh = Gridap.CellData.change_domain(u,ReferenceDomain(), Ωhact, ReferenceDomain())
  vh = Gridap.CellData.change_domain(v,ReferenceDomain(), Ωhact, ReferenceDomain())
  ∫( gvh⋅guh ) * dΩ +
  ∫( (γd/h)*vh*uh  - vh*(n_Γ⋅guh) - (n_Γ⋅gvh)*uh ) * dΓ +
  ∫( (γg*h)*jump(n_Γg⋅gvh)*jump(n_Γg⋅guh) ) * dΓg
end 

function l(v)
  gvh = Gridap.CellData.change_domain(∇(v),ReferenceDomain(), Ωhact, ReferenceDomain()) 
  vh = Gridap.CellData.change_domain(v,ReferenceDomain(), Ωhact, ReferenceDomain())
  ∫( vh*f ) * dΩ +
    ∫( (γd/h)*vh*ud - (n_Γ⋅gvh)*ud ) * dΓ
end

# FE problem
op = AffineFEOperator(a,l,U,V)
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
@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

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