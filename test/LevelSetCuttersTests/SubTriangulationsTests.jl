module SubTriangulationsTests

using Gridap
using Gridap.Geometry
using Gridap.Visualization
using GridapEmbedded.CSG
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.LevelSetCutters: initial_sub_triangulation
using GridapEmbedded.LevelSetCutters: cut_sub_triangulation


using GridapEmbedded.Interfaces

function compute_inout(a::Leaf)
  first(a.data)
end

function compute_inout(a::Node)
  cell_to_inout_1 = compute_inout(a.leftchild)
  cell_to_inout_2 = compute_inout(a.rightchild)
  op = first(a.data)
  if op  == :∪
    return _compute_inout_union.(cell_to_inout_1,cell_to_inout_2)
  elseif op == :∩
    return _compute_inout_intersection.(cell_to_inout_1,cell_to_inout_2)
  else
    @error "operation $op not implemented"
  end
end

function _compute_inout_union(inout_1,inout_2)
  if (inout_1==IN) || (inout_2==IN)
    Int8(IN)
  else
    Int8(OUT)
  end
end

function _compute_inout_intersection(inout_1,inout_2)
  if (inout_1==IN) && (inout_2==IN)
    Int8(IN)
  else
    Int8(OUT)
  end
end

R = 0.85

geo1 = disk(R,name="disk1")
geo2 = disk(R,x0=Point(1,1),name="disk2")
geo3 = disk(R,x0=Point(1,0),name="disk3")
geo5 = union(geo1,geo2)
geo4 = union(geo5,geo3)
geo6 = intersect(geo3,geo2)


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


function conversion(data)
  f,name,meta = data
  oid = objectid(f)
  ls = oid_to_ls[oid]
  cell_to_inout = ls_to_cell_to_inout[ls]
  cell_to_inout, name, meta
end

tree = replace_data(identity,conversion,get_tree(geo6))

cell_to_inout = compute_inout(tree)

celldata = ["inout_$i"=>j for (i,j) in enumerate(ls_to_cell_to_inout)]
push!(celldata,"inout"=>cell_to_inout)

write_vtk_file(UnstructuredGrid(subtrian),"subtrian2",celldata=celldata)





#st, ls_to_fst = cut_sub_triangulation(subtrian,subgeom)
#
#write_vtk_file(grid,"grid")
#writevtk(st,"st")
#
#for (i,fst) in enumerate(ls_to_fst)
#  writevtk(fst,"fst_$i")
#end

end #module
