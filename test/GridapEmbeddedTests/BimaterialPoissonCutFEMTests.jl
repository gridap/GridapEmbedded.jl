module BimaterialPoissonCutFEMTests

using Gridap
using Gridap.Arrays
import Gridap: ∇
using GridapEmbedded
using Test

# Formulation: cutfem paper section 2.1
# Stabilization coefficient like in reference [21] of this paper
# Stabilization rhs contribution associated with the flux jump at the interface also like in [21]

# Manufactured solution
const α1 = 4.0
const α2 = 3.0
u(x) = x[1] - 0.5*x[2]
∇u(x) = VectorValue( 1.0, -0.5)
q1(x) = α1*∇u(x)
q2(x) = α2*∇u(x)
j(x) = q1(x)-q2(x)
Δu(x) = 0
f1(x) = - α1*Δu(x)
f2(x) = - α2*Δu(x)
ud(x) = u(x)
∇(::typeof(u)) = ∇u

# Select geometry
const R = 0.7
geom = disc(R)
n = 30
domain = (-1,1,-1,1)
partition = (n,n)

# Setup background model
bgmodel = simplexify(CartesianDiscreteModel(domain,partition))

# Cut the background model
cutdisc = cut(bgmodel,geom)

# Setup models
model1 = DiscreteModel(cutdisc,IN)
model2 = DiscreteModel(cutdisc,OUT)

# Setup integration meshes
trian_Ω1 = Triangulation(cutdisc,IN)
trian_Ω2 = Triangulation(cutdisc,OUT)
trian_Γ = EmbeddedBoundary(cutdisc)

# Setup normal vectors
const n_Γ = get_normal_vector(trian_Γ)

# Setup cuadratures
order = 1
quad_Ω1 = CellQuadrature(trian_Ω1,2*order)
quad_Ω2 = CellQuadrature(trian_Ω2,2*order)
quad_Γ = CellQuadrature(trian_Γ,2*order)

# Setup stabilization parameters

meas_K1 = cell_measure(trian_Ω1.a,num_cells(bgmodel)) # TODO .a
meas_K2 = cell_measure(trian_Ω2.a,num_cells(bgmodel)) # TODO .a
meas_KΓ = cell_measure(trian_Γ,num_cells(bgmodel))

meas_K1_Γ = reindex(meas_K1,trian_Γ)
meas_K2_Γ = reindex(meas_K2,trian_Γ)
meas_KΓ_Γ = reindex(meas_KΓ,trian_Γ)

γ_hat = 2
const κ1 = (α2*meas_K1_Γ) ./ (α2*meas_K1_Γ .+ α1*meas_K2_Γ)
const κ2 = (α1*meas_K2_Γ) ./ (α2*meas_K1_Γ .+ α1*meas_K2_Γ)
const β = (γ_hat*meas_KΓ_Γ) ./ ( meas_K1_Γ/α1 .+ meas_K2_Γ/α2 )

# Setup FESpace

V1 = TestFESpace(
  model=model1,
  valuetype=Float64,
  reffe=:Lagrangian,
  order=order,
  conformity=:H1)

V2 = TestFESpace(
  model=model2,
  valuetype=Float64,
  reffe=:Lagrangian,
  order=order,
  conformity=:H1,
  dirichlet_tags="boundary")

U1= TrialFESpace(V1)
U2 = TrialFESpace(V2,ud)
V = MultiFieldFESpace([V1,V2])
U = MultiFieldFESpace([U1,U2])

function a_Ω1(u,v)
  u1,u2 = u
  v1,v2 = v
  α1*∇(v1)*∇(u1)
end

function a_Ω2(u,v)
  u1,u2 = u
  v1,v2 = v
  α2*∇(v2)*∇(u2)
end

function l_Ω1(v)
  v1,v2 = v
  v1*f1
end

function l_Ω2(v)
  v1,v2 = v
  v2*f2
end

#function jump(u)
#  u1,u2 = u
#  u1 - u2
#end
#
#function mean_q(u)
#  u1,u2 = u
#  κ1*α1*n_Γ*∇(u1) + κ2*α2*n_Γ*∇(u2)
#end

function a_Γ(u,v)
  # TODO: Following line not yet working
  # β*jump(v)*jump(u)  - jump(v)*mean_q(u) - mean_q(v)*jump(u)
  
  u1,u2 = u
  v1,v2 = v
  A = β*v1*u1 - β*v1*u2 - β*v2*u1 + β*v2*u2
  B = v2*κ1*α1*n_Γ*∇(u1) + v2*κ2*α2*n_Γ*∇(u2)  - v1*κ1*α1*n_Γ*∇(u1) - v1*κ2*α2*n_Γ*∇(u2)
  C = u2*κ1*α1*n_Γ*∇(v1) + u2*κ2*α2*n_Γ*∇(v2)  - u1*κ1*α1*n_Γ*∇(v1) - u1*κ2*α2*n_Γ*∇(v2)
  A + B + C
end

function l_Γ(v)
  v1,v2 = v
  v1*κ2*n_Γ*j + v2*κ1*n_Γ*j
end

# FE problem
t_Ω1 = AffineFETerm(a_Ω1,l_Ω1,trian_Ω1,quad_Ω1)
t_Ω2 = AffineFETerm(a_Ω2,l_Ω2,trian_Ω2,quad_Ω2)
t_Γ = AffineFETerm(a_Γ,l_Γ,trian_Γ,quad_Γ)
op = AffineFEOperator(U,V,t_Ω1,t_Ω2,t_Γ)
uh1, uh2 = solve(op)

# Postprocess
uh_Ω1 = restrict(uh1,trian_Ω1)
uh_Ω2 = restrict(uh2,trian_Ω2)
qh_Ω1 = α1*∇(uh_Ω1)
qh_Ω2 = α2*∇(uh_Ω2)
e1 = u - uh_Ω1
e2 = u - uh_Ω2
#writevtk(trian_Ω1,"results1",cellfields=["uh"=>uh_Ω1,"qh"=>qh_Ω1,"e"=>e1])
#writevtk(trian_Ω2,"results2",cellfields=["uh"=>uh_Ω2,"qh"=>qh_Ω2,"e"=>e2])

e1l2 = sqrt(sum(integrate(e1*e1,trian_Ω1,quad_Ω1)))
e2l2 = sqrt(sum(integrate(e2*e2,trian_Ω2,quad_Ω2)))
tol = 1.0e-9
@test e1l2 < tol
@test e2l2 < tol

end # module
