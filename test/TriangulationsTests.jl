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

cst = cut_sub_triangulation(st)
writevtk(cst,"cst")

end #module
