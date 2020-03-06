module TriangulationsTests

using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using GridapEmbedded

function doughnut(x::AbstractArray{<:Point},R,r,x0=zero(eltype(x)))
  T = eltype(eltype(x))
  npoints = length(x)
  v = zeros(T,npoints)
  for point in 1:npoints
    xi = x[point]
    vi = doughnut(xi,R,r,x0)
    v[point] = vi
  end
  v
end

@inline function doughnut(x::Point,R,r,x0)
  _x = x - x0
  (R - sqrt(_x[1]^2+_x[2]^2) )^2 + _x[3]^2 - r^2
end

function doughnut_domain(R,r)
  m = 0.1*r
  A = (R+r)+m
  B = r+m
  domain = (-A,A,-A,A,-B,B)
  domain
end

function olympic_rings(x,R,r)
  z = zero(R)
  x0 = Point(-(r+r+R),z,z)
  v1 = doughnut(x,R,r,x0)
  x0 = Point(r+r+R,z,z)
  v2 = doughnut(x,R,r,x0)
  x0 = Point(z,-R+r,z)
  v3 = doughnut(x,R,r,x0)
  x0 = Point(2*(r+r+R),-R+r,z)
  v4 = doughnut(x,R,r,x0)
  x0 = Point(3*(r+r+R),z,z)
  v5 = doughnut(x,R,r,x0)
  min.(v1,v2,v3,v4,v5)
end

function olympic_rings_domain(R,r)
  m = 0.1*r
  A = 2*(R+r)+r+m
  B = (R+r)+m
  C = r+m
  domain = (-A,A+2*(r+r+R),-B-(R-r),B,-C,C)
  domain
end

const R = 1.2
const r = 0.2

domain = olympic_rings_domain(R,r)
n = 20
partition =(8*n,4*n,n)
grid = CartesianGrid(domain,partition)
writevtk(grid,"grid")

point_to_coords = get_node_coordinates(grid)
point_to_value = olympic_rings(point_to_coords,R,r)

st = initial_sub_triangulation(grid,point_to_value)
writevtk(st,"st")

cst, ls_to_fst = cut_sub_triangulation(st)
writevtk(cst,"cst")

for (i,fst) in enumerate(ls_to_fst)
  fug = UnstructuredGrid(fst)
  quad = CellQuadrature(fug,2)
  dS = integrate(1,fug,quad)
  writevtk(fug,"fug_$i",celldata=["normals"=>fst.facet_to_normal,"dS"=>dS])
end



#domain = (0,1,0,1,0,1)
#n = 50
#partition =(n,n,n)
#grid = CartesianGrid(domain,partition)
#
#point_to_coords = get_node_coordinates(grid)
#
#const R1 = 0.7
#ls1(x) = x[1]^2 + x[2]^2 - R1
#
#const R2 = 0.7
#ls2(x) = x[1]^2 + (x[2]-1)^2 - R2
#
#const R3 = 0.7
#ls3(x) = (x[1]-1)^2 + x[2]^2 - R3
#
#point_to_ls1 = ls1.(point_to_coords)
#point_to_ls2 = ls2.(point_to_coords)
#point_to_ls3 = ls3.(point_to_coords)
#
#ls_to_point_to_value = [point_to_ls1, point_to_ls2, point_to_ls3]
#
#st = initial_sub_triangulation(grid,ls_to_point_to_value)
#writevtk(st,"st")
#
#cst, ls_to_fst = cut_sub_triangulation(st)
#writevtk(cst,"cst")
#
#for (i,fst) in enumerate(ls_to_fst)
#  fug = UnstructuredGrid(fst)
#  quad = CellQuadrature(fug,2)
#  dS = integrate(1,fug,quad)
#  writevtk(fug,"fug_$i",celldata=["normals"=>fst.facet_to_normal,"dS"=>dS])
#end



end #module
