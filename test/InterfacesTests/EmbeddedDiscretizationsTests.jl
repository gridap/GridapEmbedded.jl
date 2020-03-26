module EmbeddedDiscretizationsTests

using Gridap
using Gridap.Geometry
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters

const R = 0.7
geom = disc(R)
n = 20
partition = (n,n)

model = CartesianDiscreteModel(geom.pmin,geom.pmax,partition)

cutdisc = cut(model,geom)

model_in = DiscreteModel(cutdisc)

model_out = DiscreteModel(cutdisc,OUT)

trian_Γ = EmbeddedBoundary(cutdisc)
test_triangulation(trian_Γ)

n_Γ = get_normal_vector(trian_Γ)

V_in = TestFESpace(model=model_in,valuetype=Float64,reffe=:Lagrangian,order=1,conformity=:H1)

v_in = FEFunction(V_in,rand(num_free_dofs(V_in)))

v_in_Γ = restrict(v_in,trian_Γ)

trian = Triangulation(model)

writevtk(trian,"trian",nsubcells=10,cellfields=["v_in"=>v_in])
writevtk(trian_Γ,"trian_G",order=2,cellfields=["v_in"=>v_in_Γ,"normal"=>n_Γ])

#writevtk(cutdisc,"cutdisc")
#writevtk(model_in,"model_in")
#writevtk(model_out,"model_out")


end # module
