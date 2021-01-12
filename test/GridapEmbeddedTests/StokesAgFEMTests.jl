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

Ω_bg = Triangulation(bgmodel)
Ω = Triangulation(cutgeo)
Γd = EmbeddedBoundary(cutgeo)
Γn = BoundaryTriangulation(model,tags="boundary")

order = 2
degree = 2*order
dΩ = Measure(Ω,degree)
dΓd = Measure(Γd,degree)
dΓn = Measure(Γn,degree)

n_Γd = get_normal_vector(Γd)
n_Γn = get_normal_vector(Γn)

V_cell_fe_std = FiniteElements(PhysicalDomain(),
                               model,
                               lagrangian,
                               VectorValue{2,Float64},
                               order)
Vstd = FESpace(model,V_cell_fe_std)

V_cell_fe_ser = FiniteElements(PhysicalDomain(),
                               model,
                               lagrangian,
                               VectorValue{2,Float64},
                               order,
                               space=:S)
# RMK: we don't neet to impose continuity since
# we only use the cell dof basis / shapefuns
Vser = FESpace(model,V_cell_fe_ser,conformity=:L2)

Q_cell_fe_std = FiniteElements(PhysicalDomain(),
                               model,
                               lagrangian,
                               Float64,
                               order-1,
                               space=:P)
Qstd = FESpace(model,Q_cell_fe_std,conformity=:L2)

V = AgFEMSpace(Vstd,aggregates,Vser)
Q = AgFEMSpace(Qstd,aggregates)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

const γ = order*(order+1)

a((u,p),(v,q)) =
  ∫( ∇(v)⊙∇(u) - q*(∇⋅u) - (∇⋅v)*p ) * dΩ +
  ∫( (γ/h)*v⋅u - v⋅(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))⋅u + (p*n_Γd)⋅v + (q*n_Γd)⋅u ) * dΓd

l((v,q)) =
  ∫( v⋅f - q*g ) * dΩ +
  ∫( (γ/h)*v⋅u - (n_Γd⋅∇(v))⋅u + (q*n_Γd)⋅u ) * dΓd +
  ∫( v⋅(n_Γn⋅∇(u)) - (n_Γn⋅v)*p ) * dΓn

op = AffineFEOperator(a,l,X,Y)

uh, ph = solve(op)

eu = u - uh
ep = p - ph

l2(u) = sqrt(sum( ∫( u⊙u )*dΩ ))
h1(u) = sqrt(sum( ∫( u⊙u + ∇(u)⊙∇(u) )*dΩ ))

eu_l2 = l2(eu)
eu_h1 = h1(eu)
ep_l2 = l2(ep)

# #colors = color_aggregates(aggregates,bgmodel)
# #writevtk(Ω_bg,"trian",celldata=["aggregate"=>aggregates,"color"=>colors],cellfields=["uh"=>uh,"ph"=>ph])
# #writevtk(model,"model")
# #writevtk(Ω,"trian_O",cellfields=["uh"=>uh,"ph"=>ph])
# #writevtk(Γd,"trian_Gd",cellfields=["normal"=>n_Γd])
# #writevtk(Γn,"trian_Gn",cellfields=["normal"=>n_Γn])

tol = 1.0e-6
@test eu_l2 < tol
@test eu_h1 < tol
@test ep_l2 < tol

end # module
