module StokesCutFEMTests

using Gridap
import Gridap: ∇
using GridapEmbedded
using Test
using LinearAlgebra: tr

u(x) = VectorValue(x[1]*x[1], x[2], x[3])
∇u(x) = TensorValue(2*x[1],0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0)
Δu(x) = VectorValue(2.0,0.0,0.0)

p(x) = x[1] - x[2]
∇p(x) = VectorValue(1.0,-1.0,0.0)

f(x) = - Δu(x) + ∇p(x)
g(x) = tr(∇u(x))
ud(x) = u(x)

∇(::typeof(u)) = ∇u
∇(::typeof(p)) = ∇p

# Formulation taken from
# André Massing · Mats G. Larson · Anders Logg · Marie E. Rognes,
# A Stabilized Nitsche Fictitious Domain Method for the Stokes Problem
# J Sci Comput (2014) 61:604–628 DOI 10.1007/s10915-014-9838-9

# Select geometry
R = 0.7
L = 5.0
geo1 = tube(R,L,x0=Point(-0.5,0.0,-0.25),v=VectorValue(1.0,0.0,0.25))
box = get_metadata(geo1)

# Cut the background model
n = 8
partition = (5*n,n,2*n)
bgmodel = simplexify(CartesianDiscreteModel(box.pmin,box.pmax,partition))

# Cut the background model
cutgeo = cut(bgmodel,geo1)

# Generate the "active" model
model = DiscreteModel(cutgeo)

# Setup integration meshes
trian_Ω = Triangulation(cutgeo)
trian_Γ = EmbeddedBoundary(cutgeo)
trian_Γg = GhostSkeleton(cutgeo)
#trian_Γi = SkeletonTriangulation(model) # TODO This is not exactly the same mesh as in the paper

# Setup normal vectors
n_Γ = get_normal_vector(trian_Γ)
n_Γg = get_normal_vector(trian_Γg)

# Setup cuadratures
order = 1
quad_Ω = CellQuadrature(trian_Ω,2*order)
quad_Γ = CellQuadrature(trian_Γ,2*order)
quad_Γg = CellQuadrature(trian_Γg,2*order)
#quad_Γi = CellQuadrature(trian_Γi,2*order)

#writevtk(trian_Γi,"trian_Gi")
#kk

# Setup FESpace

V = TestFESpace(
  model=model,valuetype=VectorValue{3,Float64},reffe=:PLagrangian,
  order=order,conformity=:H1)

_Q = TestFESpace(
  model=model,valuetype=Float64,reffe=:PLagrangian,
  order=order,conformity=:H1)

# TODO
using Gridap.FESpaces: ZeroMeanFESpace
trian = Triangulation(bgmodel)
quad = CellQuadrature(trian,1)
Q = ZeroMeanFESpace(_Q,trian,quad)

U = TrialFESpace(V)
P = TrialFESpace(Q)

X = MultiFieldFESpace([U,P])
Y = MultiFieldFESpace([V,Q])

# Stabilization parameters
β0 = 0.25
β1 = 0.2
β2 = 0.1
γ = 10
h = (box.pmax-box.pmin)[1]/partition[1]

# Weak form
a_Ω(u,v) = inner(∇(u),∇(v))
b_Ω(v,p) = - (∇*v)*p
#c_Γi(p,q) = (β0*h)*jump(p)*jump(q)
c_Ω(p,q) = (β1*h^2)*∇(p)*∇(q)
a_Γ(u,v) = - (n_Γ*∇(u))*v - u*(n_Γ*∇(v)) + (γ/h)*u*v
b_Γ(v,p) = (n_Γ*v)*p
i_Γg(u,v) = (β2*h)*jump(n_Γg*∇(u))*jump(n_Γg*∇(v))
j_Γg(p,q) = (β0*h)*jump(p)*jump(q)
ϕ_Ω(q) = (β1*h^2)*∇(q)*f

function A_Ω(X,Y)
  u,p = X
  v,q = Y
  a_Ω(u,v)+b_Ω(u,q)+b_Ω(v,p)-c_Ω(p,q)
end

#function A_Γi(X,Y)
#  u,p = X
#  v,q = Y
#  -c_Γi(p,q)
#end

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
  (γ/h)*v*ud - ud*(n_Γ*∇(v)) + (q*n_Γ)*ud
  #ud*( (γ/h)*ud*v - (n_Γ*∇(v)) + q*n_Γ ) TODO
end

# FE problem

t_Ω = AffineFETerm(A_Ω,L_Ω,trian_Ω,quad_Ω)
t_Γ = AffineFETerm(A_Γ,L_Γ,trian_Γ,quad_Γ)
#t_Γi = LinearFETerm(A_Γi,trian_Γi,quad_Γi)
t_Γg = LinearFETerm(J_Γg,trian_Γg,quad_Γg)
op = AffineFEOperator(X,Y,t_Ω,t_Γ,t_Γg)#,t_Γi)
uh, ph = solve(op)

end # module
