module SubTriangulationsTests

using Gridap
using Gridap.Geometry
using Gridap.Visualization
using GridapEmbedded.CSG
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.LevelSetCutters: initial_sub_triangulation
using GridapEmbedded.LevelSetCutters: cut_sub_triangulation

R = 0.85

geo1 = disk(R,name="disk1")
geo2 = disk(R,x0=Point(1,1),name="disk2")
geo3 = disk(R,x0=Point(1,0),name="disk3")
geo4 = union(union(geo1,geo2),geo3)

n = 40
partition = (n,n)
pmin = Point(-1,-1)
pmax = Point(2,2)
grid = CartesianGrid(pmin,pmax,partition)

out = initial_sub_triangulation(grid,geo4)

subtrian, ls_to_point_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls = out

write_vtk_file(grid,"grid",celldata=[ "ls_$i"=>j for (i,j) in enumerate(ls_to_bgcell_to_inoutcut)])

write_vtk_file(UnstructuredGrid(subtrian),"subtrian",
  nodaldata=["lsv_$i"=>j for (i,j) in enumerate(ls_to_point_to_value)])

subtrian, ls_to_cell_to_inout = cut_sub_triangulation(subtrian,ls_to_point_to_value)

write_vtk_file(UnstructuredGrid(subtrian),"subtrian2",
  celldata=["inout_$i"=>j for (i,j) in enumerate(ls_to_cell_to_inout)])


#st, ls_to_fst = cut_sub_triangulation(subtrian,subgeom)
#
#write_vtk_file(grid,"grid")
#writevtk(st,"st")
#
#for (i,fst) in enumerate(ls_to_fst)
#  writevtk(fst,"fst_$i")
#end

end #module
