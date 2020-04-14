module AnalyticalGeometriesTests

using Gridap
using GridapEmbedded.CSG
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters

R = 0.7
r = 0.15
geo1 = doughnut(R,r)
test_geometry(geo1)

geo2 = doughnut(R,r,x0=Point(0.5,0.5,0.5),name="doughnut2")
test_geometry(geo2)

geo3 = union(geo1,geo2)
test_geometry(geo3)

geo4 = sphere(0.8*R)
test_geometry(geo4)

geo5 = intersect(geo3,geo4)
test_geometry(geo5)

geo6 = disk(0.8*R)
test_geometry(geo6)

geo7 = disk(0.8*R,x0=Point(1,1),name="disk2")

geo8 = union(geo6,geo7)

#using AbstractTrees
##print_tree(stdout,get_tree(geo1))
##print_tree(stdout,get_tree(geo2))
#print_tree(stdout,get_tree(geo5))
#print_tree(stdout,get_tree(geo8))

R = 1.2
r = 0.2
geo1 = olympic_rings(R,r)
box = get_metadata(geo1)
n = 7
partition =(20*n,10*n,n)
model = CartesianDiscreteModel(box.pmin,box.pmax,partition)

cutgeo = cut(model,geo1)

trian = Triangulation(model)
#writevtk(trian,"trian")

trian1 = Triangulation(cutgeo)
#writevtk(trian1,"trian1")

trian_Γ = EmbeddedBoundary(cutgeo)
#writevtk(trian_Γ,"trian_G",cellfields=["normal"=>get_normal_vector(trian_Γ)])

R = 0.7
L = 5.0
geo1 = tube(R,L,x0=Point(-0.5,0.0,-0.25),v=VectorValue(1.0,0.0,0.25))
box = get_metadata(geo1)

n = 40
partition = (n,n,n)
model = CartesianDiscreteModel(box.pmin,box.pmax,partition)

cutgeo = cut(model,geo1)

trian = Triangulation(model)
#writevtk(trian,"trian")

trian1 = Triangulation(cutgeo)
#writevtk(trian1,"trian1")

trian_Γ = EmbeddedBoundary(cutgeo)
#writevtk(trian_Γ,"trian_G",cellfields=["normal"=>get_normal_vector(trian_Γ)])

n = 40
partition = (n,n,n)
pmin = 0.8*Point(-1,-1,-1)
pmax = 0.8*Point(1,1,1)
model = CartesianDiscreteModel(pmin,pmax,partition)

R = 0.5
geo1 = cylinder(R,v=VectorValue(1,0,0))
geo2 = cylinder(R,v=VectorValue(0,1,0))
geo3 = cylinder(R,v=VectorValue(0,0,1))
geo4 = union(union(geo1,geo2),geo3)
geo5 = sphere(1)
geo6 = cube(L=1.5)
geo7 = intersect(geo6,geo5)
geo8 = setdiff(geo7,geo4)

#using AbstractTrees
#print_tree(stdout,get_tree(geo8))

cutgeo = cut(model,geo8)

trian = Triangulation(model)
#writevtk(trian,"trian")

trian1 = Triangulation(cutgeo)
#writevtk(trian1,"trian1")

trian4_Γ = EmbeddedBoundary(cutgeo,geo8,geo4)
#writevtk(trian4_Γ,"trian4_G")

trian5_Γ = EmbeddedBoundary(cutgeo,geo8,geo5)
#writevtk(trian5_Γ,"trian5_G")

trian6_Γ = EmbeddedBoundary(cutgeo,geo8,geo6)
#writevtk(trian6_Γ,"trian6_G")

end # module
