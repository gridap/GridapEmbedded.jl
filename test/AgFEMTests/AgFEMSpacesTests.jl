module AgFEMSpacesTests

using Gridap
using GridapEmbedded

const R = 1.4
geom = disk(R,x0=Point(1.0,1.7))
n = 3
partition = (n,n)

domain = (0,1,0,1)
bgmodel = CartesianDiscreteModel(domain,partition)

cutdisc = cut(bgmodel,geom)

trian = Triangulation(bgmodel)
trian_Ω = Triangulation(cutdisc)

model = DiscreteModel(cutdisc)
order = 1
V = TestFESpace(
  model=model,valuetype=Float64,reffe=:Lagrangian,
  order=order,conformity=:H1,dof_space=:physical)

vh = FEFunction(V,collect(1:num_free_dofs(V)))

cell_to_cellin = [0,0,9,8,8,9,8,8,9]

Vagg = AgFEMSpace(V,cell_to_cellin)

vhagg = FEFunction(Vagg,rand(num_free_dofs(Vagg)))

#writevtk(trian,"trian",nsubcells=10,cellfields=["vh"=>vh,"vhagg"=>vhagg],celldata=["acell"=>V.trian.oldcell_to_cell])
#writevtk(trian_Ω,"trian_0")

end # module
