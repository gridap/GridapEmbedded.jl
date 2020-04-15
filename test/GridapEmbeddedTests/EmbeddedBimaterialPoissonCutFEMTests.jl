module EmbeddedBimaterialPoissonCutFEMTests

using Gridap
using Gridap.Arrays
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
trian_Ω1 = Triangulation(cutgeo,geo1)
trian_Ω2 = Triangulation(cutgeo,geo2)
trian_Γ = EmbeddedBoundary(cutgeo,geo1,geo2)
trian_Γd = EmbeddedBoundary(cutgeo,geo3)
trian_Γg = GhostSkeleton(cutgeo,geo3)

# Setup normal vectors
n_Γ = get_normal_vector(trian_Γ)
n_Γd = get_normal_vector(trian_Γd)
n_Γg = get_normal_vector(trian_Γg)

# Setup quadratures
order = 1
quad_Ω1 = CellQuadrature(trian_Ω1,2*order)
quad_Ω2 = CellQuadrature(trian_Ω2,2*order)
quad_Γ = CellQuadrature(trian_Γ,2*order)
quad_Γd = CellQuadrature(trian_Γd,2*order)
quad_Γg = CellQuadrature(trian_Γg,2*order)

# Setup stabilization parameters

meas_K1 = cell_measure(trian_Ω1,num_cells(bgmodel))
meas_K2 = cell_measure(trian_Ω2,num_cells(bgmodel))
meas_KΓ = cell_measure(trian_Γ,num_cells(bgmodel))

meas_K1_Γ = reindex(meas_K1,trian_Γ)
meas_K2_Γ = reindex(meas_K2,trian_Γ)
meas_KΓ_Γ = reindex(meas_KΓ,trian_Γ)

γ_hat = 2
κ1 = (α2*meas_K1_Γ) ./ (α2*meas_K1_Γ .+ α1*meas_K2_Γ)
κ2 = (α1*meas_K2_Γ) ./ (α2*meas_K1_Γ .+ α1*meas_K2_Γ)
β = (γ_hat*meas_KΓ_Γ) ./ ( meas_K1_Γ/α1 .+ meas_K2_Γ/α2 )

h = (domain[2]-domain[1])/n
γd = 10.0
γg = 0.1

# Jump and mean operators for this formulation

function jump_u(u)
  u1,u2 = u
  u1 - u2
end

function mean_q(u)
  u1,u2 = u
  κ1*α1*n_Γ*∇(u1) + κ2*α2*n_Γ*∇(u2)
end

function mean_u(u)
  u1,u2 = u
  κ2*u1 + κ1*u2
end

# Setup FESpace

V1 = TestFESpace(
  model=model1, valuetype=Float64, reffe=:Lagrangian,
  order=order, conformity=:H1)

V2 = TestFESpace(
  model=model2, valuetype=Float64, reffe=:Lagrangian,
  order=order, conformity=:H1)

U1= TrialFESpace(V1)
U2 = TrialFESpace(V2)
V = MultiFieldFESpace([V1,V2])
U = MultiFieldFESpace([U1,U2])

# Weak form

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

function a_Γ(u,v)
  β*jump_u(v)*jump_u(u) - jump_u(v)*mean_q(u) - mean_q(v)*jump_u(u)
end

function a_Γd(u,v)
  u1,u2 = u
  v1,v2 = v
  (γd/h)*v2*u2 - v2*(n_Γd*∇(u2)) - (n_Γd*∇(v2))*u2
end

function a_Γg(v,u)
  u1,u2 = u
  v1,v2 = v
  (γg*h)*jump(n_Γg*∇(v2))*jump(n_Γg*∇(u2))
end

# FE problem
t_Ω1 = AffineFETerm(a_Ω1,l_Ω1,trian_Ω1,quad_Ω1)
t_Ω2 = AffineFETerm(a_Ω2,l_Ω2,trian_Ω2,quad_Ω2)
t_Γ = LinearFETerm(a_Γ,trian_Γ,quad_Γ)
t_Γd = LinearFETerm(a_Γd,trian_Γd,quad_Γd)
t_Γg = LinearFETerm(a_Γg,trian_Γg,quad_Γg)
op = AffineFEOperator(U,V,t_Ω1,t_Ω2,t_Γ,t_Γd,t_Γg)
uh1, uh2 = solve(op)

# Postprocess
uh_Ω1 = restrict(uh1,trian_Ω1)
uh_Ω2 = restrict(uh2,trian_Ω2)
qh_Ω1 = α1*∇(uh_Ω1)
qh_Ω2 = α2*∇(uh_Ω2)
writevtk(trian_Ω1,"results1",cellfields=["uh"=>uh_Ω1,"qh"=>qh_Ω1])
writevtk(trian_Ω2,"results2",cellfields=["uh"=>uh_Ω2,"qh"=>qh_Ω2])

#writevtk(model1,"model1")
#writevtk(model2,"model2")
#writevtk(trian_Ω1,"trian1")
#writevtk(trian_Ω2,"trian2")
#writevtk(trian_Γ,"trianG",cellfields=["normal"=>n_Γ])
#writevtk(trian_Γd,"trianGd",cellfields=["normal"=>n_Γd])


end # module

