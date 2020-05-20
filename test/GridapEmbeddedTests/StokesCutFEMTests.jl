module StokesCutFEMTests

using Gridap
import Gridap: ∇
using GridapEmbedded
using Test
using LinearAlgebra: tr

# Manufactured solution
u(x) = VectorValue(2*x[1],-2*x[2])
∇u(x) = TensorValue(2.0,0.0,0.0,-2.0)
Δu(x) = VectorValue(0.0,0.0)

p(x) = x[1] - x[2]
∇p(x) = VectorValue(1.0,-1.0)

#p(x) = 0
#∇p(x) = VectorValue(0,0)

∇(::typeof(u)) = ∇u
∇(::typeof(p)) = ∇p

# Forcing data
f(x) = - Δu(x) + ∇p(x)
g(x) = tr(∇u(x))
ud(x) = u(x)

# Formulation taken from
# André Massing · Mats G. Larson · Anders Logg · Marie E. Rognes,
# A Stabilized Nitsche Fictitious Domain Method for the Stokes Problem
# J Sci Comput (2014) 61:604–628 DOI 10.1007/s10915-014-9838-9

# Select geometry
R = 0.7
geo1 = disk(R)
box = get_metadata(geo1)

# Cut the background model
n = 10
partition = (n,n)
D = length(partition)
bgmodel = simplexify(CartesianDiscreteModel(box.pmin,box.pmax,partition))

# Cut the background model
cutgeo = cut(bgmodel,geo1)
cutgeo_facets = cut_facets(bgmodel,geo1)

# Generate the "active" model
model = DiscreteModel(cutgeo)

# Setup integration meshes
trian_Ω = Triangulation(cutgeo)
trian_Γ = EmbeddedBoundary(cutgeo)
trian_Γg = GhostSkeleton(cutgeo)
trian_Γi = SkeletonTriangulation(cutgeo_facets)

# Setup normal vectors
n_Γ = get_normal_vector(trian_Γ)
n_Γg = get_normal_vector(trian_Γg)
n_Γi = get_normal_vector(trian_Γi)

# Setup cuadratures
order = 1
quad_Ω = CellQuadrature(trian_Ω,2*order)
quad_Γ = CellQuadrature(trian_Γ,2*order)
quad_Γg = CellQuadrature(trian_Γg,2*order)
quad_Γi = CellQuadrature(trian_Γi,2*order)

#writevtk(trian_Γi,"trian_Gi")

# Setup FESpace

V = TestFESpace(
  model=model,valuetype=VectorValue{D,Float64},reffe=:PLagrangian,
  order=order,conformity=:H1)

Q = TestFESpace(
  model=model,valuetype=Float64,reffe=:PLagrangian,
  order=order,conformity=:H1,constraint=:zeromean)

U = TrialFESpace(V)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U,P])
Y = MultiFieldFESpace([V,Q])

# Stabilization parameters
β0 = 0.25
β1 = 0.2
β2 = 0.1
β3 = 0.05
γ = 10.0
h = (box.pmax-box.pmin)[1]/partition[1]

# Weak form
a_Ω(u,v) = inner(∇(u),∇(v))
b_Ω(v,p) = - (∇*v)*p
c_Γi(p,q) = (β0*h)*jump(p)*jump(q)
c_Ω(p,q) = (β1*h^2)*∇(p)*∇(q)
a_Γ(u,v) = - (n_Γ*∇(u))*v - u*(n_Γ*∇(v)) + (γ/h)*u*v
b_Γ(v,p) = (n_Γ*v)*p
i_Γg(u,v) = (β2*h)*jump(n_Γg*∇(u))*jump(n_Γg*∇(v))
j_Γg(p,q) = (β3*h^3)*jump(n_Γg*∇(p))*jump(n_Γg*∇(q)) + c_Γi(p,q)
ϕ_Ω(q) = (β1*h^2)*∇(q)*f

function A_Ω(X,Y)
  u,p = X
  v,q = Y
  a_Ω(u,v)+b_Ω(u,q)+b_Ω(v,p)-c_Ω(p,q)
end

function A_Γi(X,Y)
  u,p = X
  v,q = Y
  -c_Γi(p,q)
end

function A_Γ(X,Y)
  u,p = X
  v,q = Y
  a_Γ(u,v)+b_Γ(u,q)+b_Γ(v,p)
end

function J_Γg(X,Y)
  u,p = X
  v,q = Y
  i_Γg(u,v) - j_Γg(p,q) 
end

function L_Ω(Y)
  v,q = Y
  v*f - ϕ_Ω(q) - q*g
end

function L_Γ(Y)
  v,q = Y
  ud*( (γ/h)*v - n_Γ*∇(v) + q*n_Γ )
end

# FE problem
t_Ω = AffineFETerm(A_Ω,L_Ω,trian_Ω,quad_Ω)
t_Γ = AffineFETerm(A_Γ,L_Γ,trian_Γ,quad_Γ)
t_Γi = LinearFETerm(A_Γi,trian_Γi,quad_Γi)
t_Γg = LinearFETerm(J_Γg,trian_Γg,quad_Γg)
op = AffineFEOperator(X,Y,t_Ω,t_Γ,t_Γg,t_Γi)
uh, ph = solve(op)

# Postprocess
uh_Ω = restrict(uh,trian_Ω)
ph_Ω = restrict(ph,trian_Ω)
eu_Ω = u - uh_Ω
ep_Ω = p - ph_Ω
writevtk(trian_Ω,"results",
  cellfields=["uh"=>uh_Ω,"ph"=>ph_Ω,"eu"=>eu_Ω,"ep"=>ep_Ω])

# Checks
tol = 1.0e-9
l2(v) = v*v
h1(v) = inner(∇(v),∇(v))
eu_l2 = sqrt(sum(integrate(l2(eu_Ω),trian_Ω,quad_Ω)))
eu_h1 = sqrt(sum(integrate(h1(eu_Ω),trian_Ω,quad_Ω)))
ep_l2 = sqrt(sum(integrate(l2(ep_Ω),trian_Ω,quad_Ω)))
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
