module DiscreteGeometriesTests

using AbstractTrees
using Gridap
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
#print_tree(stdout,get_tree(geo4_x))
#print_tree(stdout,replace_data(d->objectid(first(d)),get_tree(geo4_x)))


end # module
