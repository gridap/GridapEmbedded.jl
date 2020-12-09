module PoissonCutFEMTests

using Gridap
using GridapEmbedded
using Test

# Manufactured solution
u(x) = x[1] + x[2] - x[3]
∇u(x) = ∇(u)(x)
f(x) = -Δ(u)(x)
ud(x) = u(x)

# Select geometry

const R = 0.7
geom = sphere(R)
n = 10
partition = (n,n,n)

##const R = 1.2
##const r = 0.4
##geom = doughnut(R,r)
##n = 20
##partition = (5*n,5*n,n)
##
##const R = 1.2
##const r = 0.2
##geom = olympic_rings(R,r)
##n = 10
##partition =(20*n,10*n,n)
##
##const R = 0.7
##const L = 5.0
##geom = tube(R,L,x0=Point(-0.5,0.0,-0.25))
##n = 30
##partition = (n,n,n)

# Setup background model
box = get_metadata(geom)
bgmodel = CartesianDiscreteModel(box.pmin,box.pmax,partition)
dp = box.pmax - box.pmin
const h = dp[1]/n

# Cut the background model
cutdisc = cut(bgmodel,geom)

# Setup integration meshes
Ω = Triangulation(cutdisc)
Γ = EmbeddedBoundary(cutdisc)
Γg = GhostSkeleton(cutdisc)

# Setup normal vectors
n_Γ = get_normal_vector(Γ)
n_Γg = get_normal_vector(Γg)

# Setup Lebesgue measures
order = 1
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
dΓg = Measure(Γg,degree)

# Setup FESpace
model = DiscreteModel(cutdisc)
V = TestFESpace(model,ReferenceFE(:Lagrangian,Float64,order),conformity=:H1)
U = TrialFESpace(V)

# Weak form Nitsche + ghost penalty (CutFEM paper Sect. 6.1)
const γd = 10.0
const γg = 0.1

a(u,v) =
  ∫( ∇(v)⋅∇(u) ) * dΩ +
  ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ +
  ∫( (γg*h)*jump(n_Γg⋅∇(v))*jump(n_Γg⋅∇(u)) ) * dΓg

l(v) =
  ∫( v*f ) * dΩ +
  ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

# FE problem
op = AffineFEOperator(a,l,U,V)
uh = solve(op)

e = u - uh

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

end # module
