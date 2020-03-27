module EmbeddedDiscretizationsTests

using Gridap
using Test
using Gridap.Geometry
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters

const R = 0.7
geom = disc(R)
n = 50
partition = (n,n)

model = CartesianDiscreteModel(geom.pmin,geom.pmax,partition)

cutdisc = cut(model,geom)

model_in = DiscreteModel(cutdisc)

model_out = DiscreteModel(cutdisc,OUT)

trian_in = Triangulation(cutdisc)
test_triangulation(trian_in)
quad_in = CellQuadrature(trian_in,1)
vol = sum(integrate(1,trian_in,quad_in))
@test abs(pi*R^2 - vol) < 1.0e-3

trian_Γ = EmbeddedBoundary(cutdisc)
test_triangulation(trian_Γ)
quad_Γ = CellQuadrature(trian_Γ,1)
surf = sum(integrate(1,trian_Γ,quad_Γ))
@test abs(surf - 2*pi*R) < 1.0e-3

trian_Γg_in = GhostSkeleton(cutdisc)

trian_Γg_out = GhostSkeleton(cutdisc,OUT)

n_Γ = get_normal_vector(trian_Γ)

V_in = TestFESpace(model=model_in,valuetype=Float64,reffe=:Lagrangian,order=1,conformity=:H1)

v_in = FEFunction(V_in,rand(num_free_dofs(V_in)))

v_in_Γ = restrict(v_in,trian_Γ)

v_in_in = restrict(v_in,trian_in)

v_in_Γg_in = restrict(v_in,trian_Γg_in)

trian = Triangulation(model)

d = mktempdir()
writevtk(trian,joinpath(d,"trian"),nsubcells=10,cellfields=["v_in"=>v_in])
writevtk(trian_Γ,joinpath(d,"trian_G"),order=2,cellfields=["v_in"=>v_in_Γ,"normal"=>n_Γ])
writevtk(trian_Γg_in,joinpath(d,"trian_Gg_in"),cellfields=["v_in"=>mean(v_in_Γg_in)])
writevtk(trian_Γg_out,joinpath(d,"trian_Gg_out"))
writevtk(trian_in,joinpath(d,"trian_in"),order=2,cellfields=["v_in"=>v_in_in])
writevtk(cutdisc,joinpath(d,"cutdisc"))
writevtk(model_in,joinpath(d,"model_in"))
writevtk(model_out,joinpath(d,"model_out"))
rm(d,recursive=true)

end # module
