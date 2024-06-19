# module PoissonAgFEMTests

using Gridap
using GridapEmbedded
using Test

u(x) = (x[1]-0.5)^2 + (x[2]-0.5)^2
f(x) = -Δ(u)(x)
ud(x) = u(x)

order = 2
n = 31
partition = (n,n)
domain = (0,1,0,1)
bgmodel = CartesianDiscreteModel(domain,partition;isperiodic=(true,true))
dp = 1
const h = dp/n

reffe = ReferenceFE(lagrangian,Float64,order)
Ω_bg = Triangulation(bgmodel)
V_bg = FESpace(Ω_bg,reffe)
φh = interpolate(x->sqrt((x[1]-0.5)^2+(x[2]-0.5)^2)-0.55,V_bg)
geo = DiscreteGeometry(φh,bgmodel)
cutgeo = cut(bgmodel,geo)

strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutgeo)

Ω_bg = Triangulation(bgmodel)
Ω_act = Triangulation(cutgeo,ACTIVE)
Ω = Triangulation(cutgeo,PHYSICAL)
Γ = EmbeddedBoundary(cutgeo)

n_Γ = get_normal_vector(Γ)

degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)

model = get_active_model(Ω_act)
Vstd = FESpace(Ω_act,FiniteElements(PhysicalDomain(),model,lagrangian,Float64,order))

V = AgFEMSpace(Vstd,aggregates)
U = TrialFESpace(V)

const γd = 10.0

a(u,v) =
  ∫( ∇(v)⋅∇(u) ) * dΩ +
  ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ

l(v) =
  ∫( v*f ) * dΩ +
  ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

op = AffineFEOperator(a,l,U,V)
uh = solve(op)

e = u - uh

l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

el2 = l2(e)
eh1 = h1(e)
ul2 = l2(uh)
uh1 = h1(uh)

colors = color_aggregates(aggregates,bgmodel)
writevtk(Ω_bg,"trian",celldata=["aggregate"=>aggregates,"color"=>colors],cellfields=["uh"=>uh])
writevtk(Ω,"trian_O",cellfields=["uh"=>uh])
writevtk(Γ,"trian_G")
@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

# end # module
