module TriangulationsTests

using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using GridapEmbedded

domain = (0,1,0,1,0,1)
n = 50
partition =(n,n,n)
grid = CartesianGrid(domain,partition)

point_to_coords = get_node_coordinates(grid)

const R1 = 0.7
ls1(x) = x[1]^2 + x[2]^2 - R1

const R2 = 0.7
ls2(x) = x[1]^2 + (x[2]-1)^2 - R2

const R3 = 0.7
ls3(x) = (x[1]-1)^2 + x[2]^2 - R3

point_to_ls1 = ls1.(point_to_coords)
point_to_ls2 = ls2.(point_to_coords)
point_to_ls3 = ls3.(point_to_coords)

ls_to_point_to_value = [point_to_ls1, point_to_ls2, point_to_ls3]

st = initial_sub_triangulation(grid,ls_to_point_to_value)
writevtk(st,"st")

cst, ls_to_fst = cut_sub_triangulation(st)
writevtk(cst,"cst")

cug = UnstructuredGrid(st)
quad = CellQuadrature(cug,2)

using LinearAlgebra
using Gridap.Arrays

cell_map = get_cell_map(cug)
q = get_coordinates(quad)
J = ∇(cell_map)
J_q = evaluate(J,q)
dV = apply(sum,apply(bcast(det),J_q))

writevtk(cug,"cug",celldata=["dV"=>dV])


fst = ls_to_fst[3]
ug1 = UnstructuredGrid(fst)
quad = CellQuadrature(ug1,2)

dS = integrate(1,ug1,quad)

writevtk(ug1,"ug1",celldata=["normals"=>fst.facet_to_normal,"dS"=>dS])


#domain = (0,1,0,1)#,0,1)
#n = 20
#partition =(n,n)#,n)
#grid = CartesianGrid(domain,partition)
#
#point_to_coords = get_node_coordinates(grid)
#
#const R1 = 1.5
#ls1(x) = x[1]^2 + x[2]^2 - R1
#
#point_to_ls1 = ls1.(point_to_coords)
#
#ls_to_point_to_value = [point_to_ls1,]
#
#st = initial_sub_triangulation(grid,ls_to_point_to_value)
#writevtk(st,"st")
#
#cst, ls_to_fst = cut_sub_triangulation(st)
#writevtk(cst,"cst")
#
#fst = ls_to_fst[1]
#ug1 = UnstructuredGrid(fst)
#writevtk(ug1,"ug1",celldata=["normals"=>fst.facet_to_normal])



end #module
