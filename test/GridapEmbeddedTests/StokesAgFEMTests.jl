module StokesAgFEMTests

using Gridap
using GridapEmbedded
using Test

u(x) = VectorValue(x[1]*x[1], x[2])
p(x) = x[1] - x[2]

f(x) = - Δ(u)(x) + ∇(p)(x)
g(x) = (∇⋅u)(x)

R = 0.3
pmin = Point(0,0)
pmax = Point(1,1)
n = 20
partition = (n,n)

geo1 = disk(R,x0=Point(0.5,0.5))
geo2 = ! geo1

bgmodel = CartesianDiscreteModel(pmin,pmax,partition)
dp = pmax - pmin
const h = dp[1]/n

cutgeo = cut(bgmodel,geo2)
model = DiscreteModel(cutgeo)

strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutgeo)

trian = Triangulation(bgmodel)
trian_Ω = Triangulation(cutgeo)
trian_Γd = EmbeddedBoundary(cutgeo)
trian_Γn = BoundaryTriangulation(model,"boundary")

order = 2
quad_Ω = CellQuadrature(trian_Ω,2*order)
quad_Γd = CellQuadrature(trian_Γd,2*order)
quad_Γn = CellQuadrature(trian_Γn,2*order)

n_Γd = get_normal_vector(trian_Γd)
n_Γn = get_normal_vector(trian_Γn)

Vstd = TestFESpace(
  reffe=:QLagrangian,
  conformity=:H1,
  valuetype=VectorValue{2,Float64},
  model=model,
  order=order,
  dof_space=:physical)

Vser = TestFESpace(
  reffe=:SLagrangian,
  conformity=:L2, # we don't neet to impose continuity since we only use the cell dof basis / shapefuns
  valuetype=VectorValue{2,Float64},
  model=model,
  order=order,
  dof_space=:physical)

Qstd = TestFESpace(
  reffe=:PLagrangian,
  conformity=:L2,
  valuetype=Float64,
  model=model,
  order=order-1,
  dof_space=:physical)

V = AgFEMSpace(Vstd,aggregates,Vser)
Q = AgFEMSpace(Qstd,aggregates)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

function A_Ω(x,y)
  u, p = x
  v, q = y
  ∇(v)⊙∇(u) - q*(∇⋅u) - (∇⋅v)*p
end

function B_Ω(y)
  v, q = y
  v⋅f - q*g
end

const γ = order*(order+1)

function A_Γd(x,y)
  u, p = x
  v, q = y
  (γ/h)*v⋅u - v⋅(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))⋅u + (p*n_Γd)⋅v + (q*n_Γd)⋅u
end

function B_Γd(y)
  v, q = y
  (γ/h)*v⋅u - (n_Γd⋅∇(v))⋅u + (q*n_Γd)⋅u
end

function B_Γn(y)
  v, q = y
  v⋅(n_Γn⋅∇(u)) - (n_Γn⋅v)*p
end

t_Ω = AffineFETerm(A_Ω,B_Ω,trian_Ω,quad_Ω)
t_Γd = AffineFETerm(A_Γd,B_Γd,trian_Γd,quad_Γd)
t_Γn = FESource(B_Γn,trian_Γn,quad_Γn)
op = AffineFEOperator(X,Y,t_Ω,t_Γd,t_Γn)
uh, ph = solve(op)

uh_Ω = restrict(uh,trian_Ω)
ph_Ω = restrict(ph,trian_Ω)

#colors = color_aggregates(aggregates,bgmodel)
#writevtk(trian,"trian",celldata=["aggregate"=>aggregates,"color"=>colors],cellfields=["uh"=>uh,"ph"=>ph])
#writevtk(model,"model")
#writevtk(trian_Ω,"trian_O",cellfields=["uh"=>uh_Ω,"ph"=>ph_Ω])
#writevtk(trian_Γd,"trian_Gd",cellfields=["normal"=>n_Γd])
#writevtk(trian_Γn,"trian_Gn",cellfields=["normal"=>n_Γn])

eu = u - uh_Ω
ep = p - ph_Ω

l2(u) = u⊙u
h1(u) = ∇(u)⊙∇(u) + l2(u)

eul2 = sqrt(sum( integrate(l2(eu),trian_Ω,quad_Ω) ))
euh1 = sqrt(sum( integrate(h1(eu),trian_Ω,quad_Ω) ))
epl2 = sqrt(sum( integrate(l2(ep),trian_Ω,quad_Ω) ))

@test eul2 < 1.e-6
@test euh1 < 1.e-6
@test epl2 < 1.e-6

end # module
