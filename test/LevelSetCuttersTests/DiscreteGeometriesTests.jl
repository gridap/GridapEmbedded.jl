module DiscreteGeometriesTests

using AbstractTrees
using Gridap
using Gridap.ReferenceFEs
using Gridap.Arrays
using GridapEmbedded.CSG
using GridapEmbedded.LevelSetCutters

R = 0.8

geo1 = disk(R,name="disk1")

geo2 = disk(R,x0=Point(1,1),name="disk2")

geo3 = union(geo1,geo2)

geo4 = intersect(geo1,geo3)

n = 10
partition = (n,n)
box1 = get_metadata(geo1)
box2 = get_metadata(geo2)
grid = CartesianGrid(box1.pmin,box2.pmax,partition)

geo4_x = discretize(geo4,grid)
test_geometry(geo4_x)

#function DiscreteGeometry test
tree1=geo4_x.tree
tree2=tree1.leftchild
level_set=tree2.data[1] #vector level set
point_to_coords=collect1d(get_node_coordinates(grid)) #model coordinates

#using new function
geo4_y=DiscreteGeometry(level_set,point_to_coords,name="")
test_geometry(geo4_y)
#print_tree(stdout,get_tree(geo4_x))
#print_tree(stdout,replace_data(d->objectid(first(d)),get_tree(geo4_x)))

#using a cellfield
domain = (0,1,0,1)
n = 10
bgmodel = CartesianDiscreteModel(domain,(n,n))

reffe = ReferenceFE(lagrangian,Float64,1)
Ω_bg = Triangulation(bgmodel)
V_bg = FESpace(Ω_bg,reffe)
φh = interpolate(x->sqrt((x[1]-0.5)^2+(x[2]-0.5)^2)-0.55,V_bg)
geo = DiscreteGeometry(φh,bgmodel)
test_geometry(geo)

end # module
