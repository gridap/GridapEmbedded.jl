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
Ω = Triangulation(cutgeo)
Γ = EmbeddedBoundary(cutgeo)
Γg = GhostSkeleton(cutgeo)
Γi = SkeletonTriangulation(cutgeo_facets)

# Setup normal vectors
n_Γ = get_normal_vector(Γ)
n_Γg = get_normal_vector(Γg)
n_Γi = get_normal_vector(Γi)

# Setup Lebesgue measures
order = 1
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
dΓg = Measure(Γg,degree)
dΓi = Measure(Γi,degree)

# Setup FESpace

reffe_u = ReferenceFE(lagrangian,VectorValue{D,Float64},order,space=:P)
reffe_p = ReferenceFE(lagrangian,Float64,order,space=:P)

V = TestFESpace(model,reffe_u,conformity=:H1)
Q = TestFESpace(model,reffe_p,conformity=:H1,constraint=:zeromean)

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
a_Ω(u,v) = ∇(u)⊙∇(v)
b_Ω(v,p) = - (∇⋅v)*p
c_Γi(p,q) = (β0*h)*jump(p)*jump(q)
c_Ω(p,q) = (β1*h^2)*∇(p)⋅∇(q)
a_Γ(u,v) = - (n_Γ⋅∇(u))⋅v - u⋅(n_Γ⋅∇(v)) + (γ/h)*u⋅v
b_Γ(v,p) = (n_Γ⋅v)*p
i_Γg(u,v) = (β2*h)*jump(n_Γg⋅∇(u))⋅jump(n_Γg⋅∇(v))
j_Γg(p,q) = (β3*h^3)*jump(n_Γg⋅∇(p))*jump(n_Γg⋅∇(q)) + c_Γi(p,q)
ϕ_Ω(q) = (β1*h^2)*∇(q)⋅f

# a((u,p),(v,q)) =
#   ∫( a_Ω(u,v)+b_Ω(u,q)+b_Ω(v,p)-c_Ω(p,q) ) * dΩ +
#   ∫( - c_Γi(p,q) ) * dΓi +
#   ∫( a_Γ(u,v)+b_Γ(u,q)+b_Γ(v,p) ) * dΓ +
#   ∫( i_Γg(u,v) - j_Γg(p,q) ) * dΓg

a((u,p),(v,q)) =
  ∫( a_Ω(u,v)+b_Ω(u,q)+b_Ω(v,p)-c_Ω(p,q) ) * dΩ +
  ∫( a_Γ(u,v)+b_Γ(u,q)+b_Γ(v,p) ) * dΓ +
  ∫( i_Γg(u,v) - j_Γg(p,q) ) * dΓg

l((v,q)) =
  ∫( v⋅f - ϕ_Ω(q) - q*g ) * dΩ +
  ∫( ud⊙( (γ/h)*v - n_Γ⋅∇(v) + q*n_Γ ) ) * dΓ

op = AffineFEOperator(a,l,X,Y)

uh, ph = solve(op)

eu = u - uh
ep = p - ph

l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))

eu_l2 = l2(eu)
eu_h1 = h1(eu)
ep_l2 = l2(ep)

# writevtk(Ω,"results",
#   cellfields=["uh"=>uh,"ph"=>ph,"eu"=>eu,"ep"=>ep])

tol = 1.0e-9
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
