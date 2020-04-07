module LevelSetCuttersTests

using Gridap
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters

R = 0.85
geo1 = disk(R,x0=Point(0,0),name="disk1")
geo2 = disk(R,x0=Point(1,1),name="disk2")
geo3 = disk(R,x0=Point(1,0),name="disk3")
geo4 = union(geo1,geo2)
geo5 = union(geo4,geo3)
geo6 = intersect(geo1,geo2)

n = 40
partition = (n,n)
pmin = Point(-1,-1)
pmax = Point(2,2)
model = CartesianDiscreteModel(pmin,pmax,partition)

cutter = LevelSetCutter()

cutgeo = cut(cutter,model,geo5)

model5 = DiscreteModel(cutgeo,geo5,(IN,CUT))
model4 = DiscreteModel(cutgeo,geo4,(IN,CUT))
model6 = DiscreteModel(cutgeo,geo6)

trian = Triangulation(model)
trian5 = Triangulation(model5)
trian4 = Triangulation(model4)
trian6 = Triangulation(model6)

writevtk(trian,"trian")
writevtk(trian5,"trian5")
writevtk(trian4,"trian4")
writevtk(trian6,"trian6")






#writevtk(cutdisc,"cutdisc")

end # module
