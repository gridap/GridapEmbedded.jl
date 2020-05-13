module PoissonAgFEMTests

using Gridap
using GridapEmbedded
using Test

u(x) = x[1] - x[2]
f(x) = -Δ(u)(x)
ud(x) = u(x)

R = 0.5
L = 0.8*(2*R)
p1 = Point(0.0,0.0)
p2 = p1 + VectorValue(L,0.0)

geo1 = disk(R,x0=p1)
geo2 = disk(R,x0=p2)
geo3 = setdiff(geo1,geo2)

a = 1.01
pmin = p1-a*R
pmax = p1+a*R

n = 30
partition = (n,n)
bgmodel = CartesianDiscreteModel(pmin,pmax,partition)
dp = pmax - pmin
const h = dp[1]/n

cutgeo = cut(bgmodel,geo3)

strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutgeo)

trian = Triangulation(bgmodel)
trian_Ω = Triangulation(cutgeo)
trian_Γ = EmbeddedBoundary(cutgeo)

n_Γ = get_normal_vector(trian_Γ)

order = 1
quad_Ω = CellQuadrature(trian_Ω,2*order)
quad_Γ = CellQuadrature(trian_Γ,2*order)

model = DiscreteModel(cutgeo)

Vstd = TestFESpace(
  model=model,valuetype=Float64,reffe=:Lagrangian,
  order=order,conformity=:H1,dof_space=:physical)

V = AgFEMSpace(Vstd,aggregates)
U = TrialFESpace(V)

const γd = 10.0
a_Ω(u,v) = ∇(v)*∇(u)
l_Ω(v) = v*f
a_Γ(u,v) = (γd/h)*v*u  - v*(n_Γ*∇(u)) - (n_Γ*∇(v))*u
l_Γ(v) = (γd/h)*v*ud - (n_Γ*∇(v))*ud

t_Ω = AffineFETerm(a_Ω,l_Ω,trian_Ω,quad_Ω)
t_Γ = AffineFETerm(a_Γ,l_Γ,trian_Γ,quad_Γ)
op = AffineFEOperator(U,V,t_Ω,t_Γ)
uh = solve(op)

uh_Ω = restrict(uh,trian_Ω)

tol = 1.0e-9
e = u - uh_Ω
el2 = sqrt(sum(integrate(e*e,trian_Ω,quad_Ω)))
@test el2 < tol
eh1 = sqrt(sum(integrate(e*e+a_Ω(e,e),trian_Ω,quad_Ω)))
@test eh1 < tol

#colors = color_aggregates(aggregates,bgmodel)
#writevtk(trian,"trian",celldata=["aggregate"=>aggregates,"color"=>colors],cellfields=["uh"=>uh])
#writevtk(trian_Ω,"trian_O",cellfields=["uh"=>uh_Ω])
#writevtk(trian_Γ,"trian_G")

end # module
