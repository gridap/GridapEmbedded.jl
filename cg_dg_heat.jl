using Gridap
using GridapEmbedded
using Gridap.ReferenceFEs
using Gridap: ∇
using Gridap.Visualization

# Setting directory for storing output files
d = "./"

# Manufactured Solution
u(x) = x[1]+x[2]
const k1 = TensorValue(1.0,0.0,0.0,0.0)
const k2 = VectorValue(0.0,1.0)
f(x) = 1.0

# Background Cartesian mesh
domain = (0,1,0,1)
n = 5
partition = (n,n)
bgmodel = CartesianDiscreteModel(domain,partition)

# Domain
const r = 0.5
geo = disk(r,x0=Point(0.5,0.5))
cutgeo = cut(bgmodel,geo)
model = DiscreteModel(cutgeo)

strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutgeo)

# Trial & Test Space
T = Float64
order = (2,1)
reffe = LagrangianRefFE(T,QUAD,order)
conf = CDConformity((CONT,DISC))
face_own_dofs = get_face_own_dofs(reffe,conf)
V0 = TestFESpace(reffe=reffe, model=model, conformity=conf,dof_space=:physical)
V = AgFEMSpace(V0,aggregates)
U = TrialFESpace(V)

# Triangulation
trian = Triangulation(bgmodel)
trian_Ω = Triangulation(cutgeo)
btrian = BoundaryTriangulation(model,[5])
strian = SkeletonTriangulation(model)
nb = get_normal_vector(btrian)
ns = get_normal_vector(strian)

# Quadrature
degree = 2 * max(order...)
quad = CellQuadrature(trian_Ω,degree)
bquad = CellQuadrature(btrian,degree)
squad = CellQuadrature(strian,degree)

# Weak formulation
a_Ω(u,v) = ∇(v)⊙(k1⋅∇(u)) + v*(k2⋅∇(u))
b_Ω(v) = v*f
t_Ω = AffineFETerm(a_Ω,b_Ω,trian_Ω,quad)

a_s(u,v) = jump(u)*(v.outward)
t_s = LinearFETerm(a_s,strian,squad)

a_d(u,v) = u*v
b_d(v) = u*v
t_d = AffineFETerm(a_d,b_d,btrian,bquad)

# Solving
op = AffineFEOperator(U,V,t_Ω,t_d,t_s)
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
el2 = sqrt(sum(integrate(l2(e),trian_Ω,quad)))
tol = 1.0e-10
@assert el2 < tol
