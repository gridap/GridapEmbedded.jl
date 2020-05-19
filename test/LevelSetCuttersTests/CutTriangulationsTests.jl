module CutTriangulationsTests

using Gridap
using Gridap.Geometry
using Gridap.Visualization
using GridapEmbedded.CSG
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.LevelSetCutters: initial_sub_triangulation
using GridapEmbedded.LevelSetCutters: cut_sub_triangulation_several_levelsets
using GridapEmbedded.LevelSetCutters: cut_sub_triangulation
using GridapEmbedded.LevelSetCutters: cut_sub_triangulation_with_boundary

R = 0.85

geo1 = disk(R,name="disk1")
geo2 = disk(R,x0=Point(1,1),name="disk2")
geo3 = disk(R,x0=Point(1,0),name="disk3")
geo5 = union(geo1,geo2)
geo4 = union(geo5,geo3)

n = 40
partition = (n,n)
pmin = Point(-1,-1)
pmax = Point(2,2)
grid = CartesianGrid(pmin,pmax,partition)

out = initial_sub_triangulation(grid,geo4)

subtrian0, ls_to_point_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls = out

point_to_value = first(ls_to_point_to_value)
subtrian, cell_to_inout = cut_sub_triangulation(subtrian0,point_to_value)

celldata1 = ["inout"=>cell_to_inout]
celldata2 = ["bgcell"=>subtrian.cell_to_bgcell]
celldata = vcat(celldata1,celldata2)

writevtk(UnstructuredGrid(subtrian),"subtrian",celldata=celldata)

subtrian, cell_to_inout, subtrian_b = cut_sub_triangulation_with_boundary(subtrian0,point_to_value)

celldata1 = ["inout"=>cell_to_inout]
celldata2 = ["bgcell"=>subtrian.cell_to_bgcell]
celldata = vcat(celldata1,celldata2)

writevtk(UnstructuredGrid(subtrian),"subtrian1",celldata=celldata)
writevtk(subtrian_b,"subtrian1_b")

subtrian, ls_to_cell_to_inoutcut = cut_sub_triangulation_several_levelsets(subtrian0,ls_to_point_to_value)

celldata1 = ["inout_$i"=>cell_to_inout for (i,cell_to_inout) in enumerate(ls_to_cell_to_inoutcut)]
celldata2 = ["bgcell"=>subtrian.cell_to_bgcell]
celldata = vcat(celldata1,celldata2)

writevtk(UnstructuredGrid(subtrian),"subtrian2",celldata=celldata)

end # module
