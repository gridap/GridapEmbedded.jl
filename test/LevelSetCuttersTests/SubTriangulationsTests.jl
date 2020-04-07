module SubTriangulationsTests

using Gridap
using Gridap.Geometry
using Gridap.Visualization
using GridapEmbedded.CSG
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.LevelSetCutters: initial_sub_triangulation
using GridapEmbedded.LevelSetCutters: cut_sub_triangulation
using GridapEmbedded.Interfaces: merge_facet_sub_triangulations

R = 0.85

geo1 = disk(R,name="disk1")
geo2 = disk(R,x0=Point(1,1),name="disk2")
geo3 = disk(R,x0=Point(1,0),name="disk3")
geo5 = union(geo1,geo2)
geo4 = union(geo5,geo3)
geo6 = intersect(geo3,geo2)
geo7 = setdiff(geo3,geo2)


n = 40
partition = (n,n)
pmin = Point(-1,-1)
pmax = Point(2,2)
grid = CartesianGrid(pmin,pmax,partition)

out = initial_sub_triangulation(grid,geo4)

subtrian, ls_to_point_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls = out

#write_vtk_file(grid,"grid",celldata=[ "ls_$i"=>j for (i,j) in enumerate(ls_to_bgcell_to_inoutcut)])
#
#write_vtk_file(UnstructuredGrid(subtrian),"subtrian",
#  nodaldata=["lsv_$i"=>j for (i,j) in enumerate(ls_to_point_to_value)])

subtrian, ls_to_cell_to_inout, fst, ls_to_facet_to_inout = cut_sub_triangulation(subtrian,ls_to_point_to_value)

celldata = ["inout_$i"=>j for (i,j) in enumerate(ls_to_cell_to_inout)]

#write_vtk_file(UnstructuredGrid(subtrian),"subtrian2",celldata=celldata)



tree = get_tree(geo7)
function conversion(data)
  f,name,meta = data
  oid = objectid(f)
  ls = oid_to_ls[oid]
  cell_to_inoutcut = ls_to_facet_to_inout[ls]
  cell_to_inoutcut, name, meta
end
newtree = replace_data(identity,conversion,tree)

using GridapEmbedded.Interfaces: compute_inoutcut

facet_to_inout = compute_inoutcut(newtree)

celldata1 = ["inout_$k"=>j for (k,j) in enumerate(ls_to_facet_to_inout)]
celldata2 = ["normal"=> fst.facet_to_normal,"inout"=>facet_to_inout]
celldata = vcat(celldata1,celldata2)
write_vtk_file(UnstructuredGrid(fst),"fst", celldata=celldata)


#function compute_inoutboundary(a::Leaf)
#  first(a.data)
#end
#
#function compute_inoutboundary(a::Node)
#  cell_to_inoutcut_1 = compute_inoutboundary(a.leftchild)
#  cell_to_inoutcut_2 = compute_inoutboundary(a.rightchild)
#  op = first(a.data)
#  if op  == :∪
#    return _compute_inoutboundary_union.(cell_to_inoutcut_1,cell_to_inoutcut_2)
#  elseif op == :∩
#    return _compute_inoutboundary_intersection.(cell_to_inoutcut_1,cell_to_inoutcut_2)
#  elseif op == :-
#    return _compute_inoutboundary_setdiff.(cell_to_inoutcut_1,cell_to_inoutcut_2)
#  else
#    @error "operation $op not implemented"
#  end
#end








end #module
