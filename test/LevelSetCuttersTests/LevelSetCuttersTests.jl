module LevelSetCuttersTests

using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: compute_bgcell_to_inoutcut
using GridapEmbedded.LevelSetCutters

R = 0.85
geo1 = disk(R,x0=Point(0,0),name="disk1")
geo2 = disk(R,x0=Point(1,1),name="disk2")
geo3 = disk(R,x0=Point(1,0),name="disk3")
geo4 = union(geo1,geo2)
geo5 = union(geo4,geo3)
geo6 = intersect(geo1,geo2)
geo7 = intersect(geo4,geo3)
geo8 = setdiff(geo1,geo2)

n = 40
partition = (n,n)
pmin = Point(-1,-1)
pmax = Point(2,2)
model = CartesianDiscreteModel(pmin,pmax,partition)
trian = Triangulation(model)
trian_facets = Triangulation(ReferenceFE{1},model)

bgcell_to_inoutcut = compute_bgcell_to_inoutcut(model,geo4)
#writevtk(trian,"trian",celldata=["inoutcut"=>bgcell_to_inoutcut])

bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(model,geo4)
#writevtk(trian_facets,"trian_facets",celldata=["inoutcut"=>bgfacet_to_inoutcut])

cutter = LevelSetCutter()

cutgeo = cut(cutter,model,geo5)

model5 = DiscreteModel(cutgeo,geo5,(IN,CUT))
model4 = DiscreteModel(cutgeo,geo4,(IN,CUT))
model6 = DiscreteModel(cutgeo,geo6)

bgcell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geo4)

trian5 = Triangulation(model5)
trian4 = Triangulation(model4)
trian6 = Triangulation(model6)

#writevtk(trian,"trian",celldata=["inoutcut"=>bgcell_to_inoutcut])
#writevtk(trian5,"trian5")
#writevtk(trian4,"trian4")
#writevtk(trian6,"trian6")

trian4cut = Triangulation(cutgeo,geo4,CUTIN)
#writevtk(trian4cut,"trian4cut")

trian4cutin = Triangulation(cutgeo,geo4,(CUTIN,IN))
#writevtk(trian4cutin,"trian4cutin")

trian6cutin = Triangulation(cutgeo,geo6,(CUTIN,IN))
#writevtk(trian6cutin,"trian6cutin")

trian7cutin = Triangulation(cutgeo,geo7,(CUTIN,IN))
#writevtk(trian7cutin,"trian7cutin")

trian8cutin = Triangulation(cutgeo,geo8,(CUTIN,IN))
#writevtk(trian8cutin,"trian8cutin")

trian5_Γ = EmbeddedBoundary(cutgeo,geo5)
#writevtk(trian5_Γ,"trian5_G")

trian8_Γ = EmbeddedBoundary(cutgeo,geo8)
#writevtk(trian8_Γ,"trian8_G")

trian82_Γ = EmbeddedBoundary(cutgeo,geo8,geo2)
#writevtk(trian82_Γ,"trian82_G")


#writevtk(cutdisc,"cutdisc")


n = 20
partition = (n,n)
pmin = Point(0,0)
pmax = Point(1,1)
model = CartesianDiscreteModel(pmin,pmax,partition)

R = 0.44
geo1 = disk(R,x0=Point(0.0,0.5),name="disk1")
geo2 = disk(R,x0=Point(0.5,1.0),name="disk2")
geo3 = union(geo1,geo2,name="domain")

cutgeo = cut(cutter,model,geo3)
cutgeo_facets = cut_facets(cutter,model,geo3)

writevtk(cutgeo_facets.subfacets,"subfacets")

trian = Triangulation(model)
trian3 = Triangulation(cutgeo,geo3)
btrian3in = BoundaryTriangulation(cutgeo_facets,"boundary",geo3,IN)
btrian3 = BoundaryTriangulation(cutgeo_facets,"boundary",geo3,(CUTIN,IN))

writevtk(trian,"trian")
writevtk(trian3,"trian3")
writevtk(btrian3in,"btrian3in",
  celldata=["bgcell"=>get_cell_id(btrian3in)],
  cellfields=["normal"=>get_normal_vector(btrian3in)])
writevtk(btrian3,"btrian3",
  celldata=["bgcell"=>get_cell_id(btrian3)],
  cellfields=["normal"=>get_normal_vector(btrian3)])

#
#btrian

end # module
