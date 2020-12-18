module BimaterialPoissonCutFEMTests

using Gridap
import Gridap: ∇
using GridapEmbedded
using Test

# Formulation: cutfem paper section 2.1
# Stabilization coefficient like in reference [21] of this paper
# rhs contribution associated with the flux jump at the interface also like in [21]

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
geo1 = disk(R)
geo2 = ! geo1
n = 30
domain = (-1,1,-1,1)
partition = (n,n)

# Setup background model
bgmodel = simplexify(CartesianDiscreteModel(domain,partition))

# Cut the background model
cutgeo = cut(bgmodel,union(geo1,geo2))

# Setup models
model1 = DiscreteModel(cutgeo,geo1)
model2 = DiscreteModel(cutgeo,geo2)

# Setup integration meshes
Ω1 = Triangulation(cutgeo,geo1)
Ω2 = Triangulation(cutgeo,geo2)
Γ = EmbeddedBoundary(cutgeo,geo1,geo2)

# Setup normal vectors
const n_Γ = get_normal_vector(Γ)

# Setup Lebesgue measures
order = 1
degree = 2*order
dΩ1 = Measure(Ω1,degree)
dΩ2 = Measure(Ω2,degree)
dΓ = Measure(Γ,degree)

# Setup FESpace

V1 = TestFESpace(model1,ReferenceFE(lagrangian,Float64,order),conformity=:H1)

V2 = TestFESpace(model2,
                 ReferenceFE(lagrangian,Float64,order),
                 conformity=:H1,
                 dirichlet_tags="boundary")

U1 = TrialFESpace(V1)
U2 = TrialFESpace(V2,ud)

V = MultiFieldFESpace([V1,V2])
U = MultiFieldFESpace([U1,U2])

# Setup stabilization parameters

meas_K1 = get_cell_measure(Ω1)
meas_K2 = get_cell_measure(Ω2)
meas_KΓ = get_cell_measure(Γ)

# meas_K1_Γ = lazy_map(Reindex(meas_K1),get_cell_to_bgcell(Γ))
# meas_K2_Γ = lazy_map(Reindex(meas_K2),get_cell_to_bgcell(Γ))
# meas_KΓ_Γ = lazy_map(Reindex(meas_KΓ),get_cell_to_bgcell(Γ))

#writevtk(model1,"model1")
#writevtk(model2,"model2")
#writevtk(Ω1,"trian1")
#writevtk(Ω2,"trian2")
#writevtk(Γ,"trianG",
#  celldata=["K1"=>meas_K1_Γ,"K2"=>meas_K2_Γ,"KG"=>meas_KΓ_Γ],
#  cellfields=["normal"=>n_Γ])

γ_hat = 2
const κ1 = (α2*meas_K1) ./ (α2*meas_K1 .+ α1*meas_K2)
const κ2 = (α1*meas_K2) ./ (α2*meas_K1 .+ α1*meas_K2)
const β = (γ_hat*meas_KΓ) ./ ( meas_K1/α1 .+ meas_K2/α2 )

# Jump and mean operators for this formulation

jump_u(u1,u2) = u1 - u2
mean_q(u1,u2) = κ1*α1*∇(u1) + κ2*α2*∇(u2)
mean_u(u1,u2) = κ2*u1 + κ1*u2

# Weak form

a((u1,u2),(v1,v2)) =
  ∫( α1*∇(v1)⋅∇(u1) ) * dΩ1 + ∫( α2*∇(v2)⋅∇(u2) ) * dΩ2 +
  ∫( β*jump_u(v1,v2)*jump_u(u1,u2)
     - n_Γ⋅mean_q(u1,u2)*jump_u(v1,v2)
     - n_Γ⋅mean_q(v1,v2)*jump_u(u1,u2) ) * dΓ

l((v1,v2)) =
  ∫( v1*f1 ) * dΩ1 + ∫( v2*f2 ) * dΩ2 +
  ∫( mean_u(v1,v2)*(n_Γ⋅j) ) * dΓ

op = AffineFEOperator(a,l,U,V)
uh1, uh2 = solve(op)
uh = (uh1,uh2)

e1 = u - uh1
e2 = u - uh2
e = (e1,e2)

l2((u1,u2)) = sqrt( sum( ∫( u1*u1 )*dΩ1 ) + sum( ∫( u2*u2 )*dΩ2 ) )
h1((u1,u2)) = sqrt( sum( ∫( u1*u1 + ∇(u1)⋅∇(u1) )*dΩ1 ) +
                    sum( ∫( u2*u2 + ∇(u2)⋅∇(u2) )*dΩ2 ) )

el2 = l2(e)
eh1 = h1(e)
ul2 = l2(uh)
uh1 = h1(uh)

# qh1 = α1*∇(uh1)
# qh2 = α2*∇(uh2)
# writevtk(Ω1,"results1",cellfields=["uh"=>uh1,"qh"=>qh1,"e"=>e1])
# writevtk(Ω2,"results2",cellfields=["uh"=>uh2,"qh"=>qh2,"e"=>e2])
@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

end # module
