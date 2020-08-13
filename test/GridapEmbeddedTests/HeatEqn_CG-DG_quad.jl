using Gridap
using GridapEmbedded
using Gridap.ReferenceFEs
using Gridap: ∇
using Gridap.Visualization
using GridapEmbedded.LevelSetCutters

# Setting directory for storing output files
d = "./"

# Manufactured Solution
u(x) = x[1]+x[2]
const k1 = TensorValue(1.0,0.0,0.0,0.0)
const k2 = VectorValue(0.0,1.0)
f(x) = 1.0

# Background Cartesian mesh
domain = (-1,1,-2,2)
n_x = 5
n_t = 5
partition = (n_x,n_t)
bgmodel = CartesianDiscreteModel(domain,partition)

# Domain
geo = par(L=1,M=1)
cutgeo = cut(bgmodel,geo)
model = DiscreteModel(cutgeo)

# AgFEM strategy
strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutgeo)

# Trial & Test Space
T = Float64
order = (2,1)
reffe = LagrangianRefFE(T,QUAD,order)
conf = CDConformity((CONT,DISC))
V0 = TestFESpace(reffe=reffe, model=model, conformity=conf,dof_space=:physical)
V = AgFEMSpace(V0,aggregates)
U = TrialFESpace(V)

# Triangulation
trian_Ω = Triangulation(cutgeo)
trian_Γ = EmbeddedBoundary(cutgeo)
strian = GhostSkeleton(cutgeo)
n_Γ = get_normal_vector(trian_Γ)
ns = get_normal_vector(strian)
fi = joinpath(d,"strian")
writevtk(strian,fi,cellfields=["normal"=>ns])

# Quadrature
degree = 2 * max(order...)
quad_Ω = CellQuadrature(trian_Ω,degree)
quad_Γ = CellQuadrature(trian_Γ,degree)
squad = CellQuadrature(strian,degree)

# Weak formulation
a_Ω(u,v) = ∇(v)⊙(k1⋅∇(u)) + v*(k2⋅∇(u))
b_Ω(v) = v*f
t_Ω = AffineFETerm(a_Ω,b_Ω,trian_Ω,quad_Ω)

a_s(u,v) = jump(u)*(v.outward)
t_s = LinearFETerm(a_s,strian,squad)

const γd = 10.0
h = 1/n_x
n_γ = k1⋅n_Γ
a_Γ(u,v) = (γd/h)*v*u  - v*(n_γ⋅(k1⋅∇(u))) - (n_γ⋅(k1⋅∇(v)))*u
l_Γ(v) = (γd/h)*v*u - (n_γ⋅(k1⋅∇(v)))*u
t_Γ = AffineFETerm(a_Γ,l_Γ,trian_Γ,quad_Γ)

# Solving
op = AffineFEOperator(U,V,t_Ω,t_Γ,t_s)
uh = solve(op)
uh_Ω = restrict(uh,trian_Ω)
fi = joinpath(d,"results")
writevtk(trian_Ω,fi,cellfields=["uh"=>uh_Ω])

# Error
e = u-uh_Ω
fi = joinpath(d,"error")
writevtk(trian_Ω,fi,cellfields=["e" => e])

# L2 Error
l2(u) = u*u
el2 = sqrt(sum(integrate(l2(e),trian_Ω,quad_Ω)))
tol = 1.0e-10
@assert el2 < tol

trian = Triangulation(bgmodel)
colors = color_aggregates(aggregates,bgmodel)
writevtk(trian,"trian",celldata=["aggregate"=>aggregates,"color"=>colors],cellfields=["uh"=>uh])
