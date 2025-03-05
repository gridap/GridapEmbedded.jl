module TraceFEMTests

using Gridap
using GridapEmbedded
using GridapEmbedded.Interfaces
using Test

# Select geometry
const R = 0.7
geom = disk(R)
n = 10
partition = (n,n)
ud(x) = x[1] - x[2]
f(x) = ud(x)

# Setup background model
box = get_metadata(geom)
dp = box.pmax - box.pmin
const h = dp[1]/n
bgmodel = simplexify(CartesianDiscreteModel(box.pmin,box.pmax,partition))
Ω_bg = Triangulation(bgmodel)

# Cut Geometry
cutgeom = cut(bgmodel,geom)
Ωc = Triangulation(cutgeom,CUT,geom)
Γ = EmbeddedBoundary(cutgeom,geom)
Γg = GhostSkeleton(cutgeom,CUT,geom)

order=1
V = TestFESpace(Ωc,ReferenceFE(lagrangian,Float64,order),conformity=:H1)
U = TrialFESpace(V)

n_Γ = get_normal_vector(Γ)
n_Γg = get_normal_vector(Γg)

dΩc = Measure(Ωc,2)
dΓ = Measure(Γ,2)
dΓg = Measure(Γg,2)

γd = 10
γg = 0.1

a(u,v)= ∫( (γd/h)*v*u   ) * dΓ +
        ∫( (γg*h)*jump(n_Γg⋅∇(v))*jump(n_Γg⋅∇(u))) * dΓg

b(v) =  ∫( (γd/h)*v*ud  ) * dΓ

op = AffineFEOperator(a,b,U,V)

uh = solve(op)

e = uh-ud
l2(v) = √(∑(∫(v*v)dΓ))

tol = 1.0e-10
@test l2(e) < tol

#writevtk(Γ,"u_Γ",cellfields=["uh"=>uh])

end
