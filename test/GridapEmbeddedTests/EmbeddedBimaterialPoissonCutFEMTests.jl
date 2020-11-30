module EmbeddedBimaterialPoissonCutFEMTests

using Gridap
using GridapEmbedded

# Material Parameters and loads
α1 = 1.0
α2 = 5.0
f1(x) = 10
f2(x) = 10

# Select geometry
R = 0.7
geo1 = disk(0.4*R,x0=Point(0.5,0.1))
geo4 = disk(R,x0=Point(0.2,0.0))
geo5 = disk(0.8*R,x0=-Point(0.3,0.0))
geo3 = union(geo4,geo5)
geo2 = setdiff(geo3,geo1)

n = 60
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
Γd = EmbeddedBoundary(cutgeo,geo3)
Γg = GhostSkeleton(cutgeo,geo3)

# Setup normal vectors
n_Γ = get_normal_vector(Γ)
n_Γd = get_normal_vector(Γd)
n_Γg = get_normal_vector(Γg)

# Setup quadratures
order = 1
degree = 2*order
dΩ1 = LebesgueMeasure(Ω1,degree)
dΩ2 = LebesgueMeasure(Ω2,degree)
dΓ = LebesgueMeasure(Γ,degree)
dΓd = LebesgueMeasure(Γd,degree)
dΓg = LebesgueMeasure(Γg,degree)

# Setup FESpace

V1 = TestFESpace(model1,ReferenceFE(:Lagrangian,Float64,order),conformity=:H1)
V2 = TestFESpace(model2,ReferenceFE(:Lagrangian,Float64,order),conformity=:H1)

U1 = TrialFESpace(V1)
U2 = TrialFESpace(V2)
V = MultiFieldFESpace([V1,V2])
U = MultiFieldFESpace([U1,U2])

# Setup stabilization parameters

meas_K1 = get_cell_measure(Ω1)
meas_K2 = get_cell_measure(Ω2)
meas_KΓ = get_cell_measure(Γ)

γ_hat = 2
κ1 = (α2*meas_K1) ./ (α2*meas_K1 .+ α1*meas_K2)
κ2 = (α1*meas_K2) ./ (α2*meas_K1 .+ α1*meas_K2)
β = (γ_hat*meas_KΓ) ./ ( meas_K1/α1 .+ meas_K2/α2 )

h = (domain[2]-domain[1])/n
γd = 10.0
γg = 0.1

# Jump and mean operators for this formulation

jump_u(u1,u2) = u1 - u2
mean_q(u1,u2) = κ1*α1*∇(u1) + κ2*α2*∇(u2)
mean_u(u1,u2) = κ2*u1 + κ1*u2

# Weak form

a((u1,u2),(v1,v2)) =
  ∫( α1*∇(v1)⋅∇(u1) ) * dΩ1 + ∫( α2*∇(v2)⋅∇(u2) ) * dΩ2 +
  ∫( β*jump_u(v1,v2)*jump_u(u1,u2)
     - n_Γ⋅mean_q(u1,u2)*jump_u(v1,v2)
     - n_Γ⋅mean_q(v1,v2)*jump_u(u1,u2) ) * dΓ +
  ∫( (γd/h)*v2*u2 - v2*(n_Γd⋅∇(u2)) - (n_Γd⋅∇(v2))*u2 ) * dΓd +
  ∫( (γg*h)*jump(n_Γg⋅∇(v2))*jump(n_Γg⋅∇(u2)) ) * dΓg

l((v1,v2)) =
  ∫( v1*f1 ) * dΩ1 + ∫( v2*f2 ) * dΩ2

# FE problem
op = AffineFEOperator(a,l,U,V)
uh1, uh2 = solve(op)
uh = (uh1,uh2)

# Postprocess
qh1 = α1*∇(uh1)
qh2 = α2*∇(uh2)
writevtk(Ω1,"results1",cellfields=["uh"=>uh1,"qh"=>qh1])
writevtk(Ω2,"results2",cellfields=["uh"=>uh2,"qh"=>qh2])

#writevtk(model1,"model1")
#writevtk(model2,"model2")
#writevtk(Ω1,"trian1")
#writevtk(Ω2,"trian2")
#writevtk(Γ,"trianG",cellfields=["normal"=>n_Γ])
#writevtk(Γd,"trianGd",cellfields=["normal"=>n_Γd])

end # module
