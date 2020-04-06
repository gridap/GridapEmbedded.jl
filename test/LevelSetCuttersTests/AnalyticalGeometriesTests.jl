module AnalyticalGeometriesTests

using Gridap
using GridapEmbedded.CSG
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

end # module
