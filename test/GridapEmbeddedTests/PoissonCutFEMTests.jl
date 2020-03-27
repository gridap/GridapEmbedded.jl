module PoissonCutFEMTests

using Gridap
import Gridap: ∇
using GridapEmbedded
using Test

# Manufactured solution
u(x) = x[1] + x[2] - x[3]
∇u(x) = VectorValue( 1, 1, -1)
Δu(x) = 0
f(x) = - Δu(x)
ud(x) = u(x)
∇(::typeof(u)) = ∇u

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
bgmodel = CartesianDiscreteModel(geom.pmin,geom.pmax,partition)
dp = geom.pmax - geom.pmin
const h = dp[1]/n

# Cut the background model
cutdisc = cut(bgmodel,geom)

# Setup integration meshes
trian_Ω = Triangulation(cutdisc)
trian_Γ = EmbeddedBoundary(cutdisc)
trian_Γg = GhostSkeleton(cutdisc)

# Setup normal vectors
n_Γ = get_normal_vector(trian_Γ)
n_Γg = get_normal_vector(trian_Γg)

# Setup cuadratures
order = 1
quad_Ω = CellQuadrature(trian_Ω,2*order)
quad_Γ = CellQuadrature(trian_Γ,2*order)
quad_Γg = CellQuadrature(trian_Γg,2*order)

# Setup FESpace
model = DiscreteModel(cutdisc)
V = TestFESpace(model=model,valuetype=Float64,reffe=:Lagrangian,order=order,conformity=:H1)
U = TrialFESpace(V)

# Weak form Nitsche + ghost penalty (CutFEM paper Sect. 6.1)
const γd = 10.0
const γg = 0.1
a_Ω(u,v) = ∇(v)*∇(u)
l_Ω(v) = v*f
a_Γ(u,v) = (γd/h)*v*u  - v*(n_Γ*∇(u)) - (n_Γ*∇(v))*u
l_Γ(v) = (γd/h)*v*ud - (n_Γ*∇(v))*ud
a_Γg(v,u) = (γg*h)*jump(n_Γg*∇(v))*jump(n_Γg*∇(u))

# FE problem
t_Ω = AffineFETerm(a_Ω,l_Ω,trian_Ω,quad_Ω)
t_Γ = AffineFETerm(a_Γ,l_Γ,trian_Γ,quad_Γ)
t_Γg = LinearFETerm(a_Γg,trian_Γg,quad_Γg)
op = AffineFEOperator(U,V,t_Ω,t_Γ,t_Γg)
uh = solve(op)

# Postprocess
uh_Ω = restrict(uh,trian_Ω)
#writevtk(trian_Ω,"results",cellfields=["uh"=>uh_Ω])

tol = 1.0e-9
e = u - uh_Ω
el2 = sqrt(sum(integrate(e*e,trian_Ω,quad_Ω)))
@test el2 < tol
##eh1 = sqrt(sum(integrate(e*e+a_Ω(e,e),trian_Ω,quad_Ω)))
##@test eh1 < tol

end # module
