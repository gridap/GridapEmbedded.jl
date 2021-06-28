module AnalyticalGeometriesTests

using Gridap
using Gridap.Geometry
using GridapEmbedded.CSG
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters
using Test

function tmpdir(f::Function)
 d = mktempdir()
 try 
   f(d)
 finally
   rm(d,recursive=true)
 end
end

tmpdir() do d

geo0 = AnalyticalGeometry(x->.9*x[1]^2+1.1*x[2]^2-0.6)
model = CartesianDiscreteModel((-1,1,-1,1),(10,10))
cutgeo = cut(model,geo0)
Γ = EmbeddedBoundary(cutgeo)
writevtk(Γ,joinpath(d,"user"))

geo0 = popcorn()
test_geometry(geo0)

box = get_metadata(geo0)
partition = (40,40,40)
model = CartesianDiscreteModel(box.pmin,box.pmax,partition)
cutgeo = cut(model,geo0)
Γ = EmbeddedBoundary(cutgeo)
writevtk(Γ,joinpath(d,"popcorn"))

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

geo1 = disk(0.8*R)
geo2 = ! geo1
box = get_metadata(geo1)
partition = (10,10)
model = CartesianDiscreteModel(box.pmin,box.pmax,partition)

cutgeo = cut(model,geo2)

trian = Triangulation(model)
writevtk(trian,joinpath(d,"trian"))

trian1 = Triangulation(cutgeo)
test_triangulation(trian1)
writevtk(trian1,joinpath(d,"trian1"))

trian_Γ = EmbeddedBoundary(cutgeo)
test_triangulation(trian_Γ)
writevtk(trian_Γ,joinpath(d,"trian_G"),cellfields=["normal"=>get_normal_vector(trian_Γ)])

geo2 = disk(0.8*R,name="geo2")
box = get_metadata(geo2)
geo3 = disk(0.4*R,x0=Point(0.8*R,0.0),name="geo3")
geo1 = setdiff(geo2,geo3,name="geo1")

partition = (10,10)
model = CartesianDiscreteModel(box.pmin,box.pmax,partition)

cutgeo = cut(model,geo1)

trian = Triangulation(model)
writevtk(trian,joinpath(d,"trian"))

trian1 = Triangulation(cutgeo,"geo1")
writevtk(trian1,joinpath(d,"trian1"))

trian_Γ = EmbeddedBoundary(cutgeo,"geo1")
writevtk(trian_Γ,joinpath(d,"trian_G"),cellfields=["normal"=>get_normal_vector(trian_Γ)])

trian_Γ13 = EmbeddedBoundary(cutgeo,"geo1","geo3")
writevtk(trian_Γ13,joinpath(d,"trian_G13"),cellfields=["normal"=>get_normal_vector(trian_Γ13)])

R = 1.2
r = 0.2
geo1 = olympic_rings(R,r)
box = get_metadata(geo1)
n = 7
partition =(20*n,10*n,n)
model = CartesianDiscreteModel(box.pmin,box.pmax,partition)

cutgeo = cut(model,geo1)

trian = Triangulation(model)
writevtk(trian,joinpath(d,"trian"))

trian1 = Triangulation(cutgeo)
writevtk(trian1,joinpath(d,"trian1"))

trian_Γ = EmbeddedBoundary(cutgeo)
writevtk(trian_Γ,joinpath(d,"trian_G"),cellfields=["normal"=>get_normal_vector(trian_Γ)])

R = 0.7
L = 5.0
geo1 = tube(R,L,x0=Point(-0.5,0.0,-0.25),v=VectorValue(1.0,0.0,0.25))
box = get_metadata(geo1)

n = 40
partition = (n,n,n)
model = CartesianDiscreteModel(box.pmin,box.pmax,partition)

cutgeo = cut(model,geo1)

trian = Triangulation(model)
writevtk(trian,joinpath(d,"trian"))

trian1 = Triangulation(cutgeo)
writevtk(trian1,joinpath(d,"trian1"))

trian_Γ = EmbeddedBoundary(cutgeo)
writevtk(trian_Γ,joinpath(d,"trian_G"),cellfields=["normal"=>get_normal_vector(trian_Γ)])

n = 40
partition = (n,n,n)
pmin = 0.8*Point(-1,-1,-1)
pmax = 0.8*Point(1,1,1)
model = CartesianDiscreteModel(pmin,pmax,partition)

R = 0.5
geo1 = cylinder(R,v=VectorValue(1,0,0))
geo2 = cylinder(R,v=VectorValue(0,1,0))
geo3 = cylinder(R,v=VectorValue(0,0,1))
geo4 = union(union(geo1,geo2),geo3)
geo5 = sphere(1)
geo6 = cube(L=1.5)
geo7 = intersect(geo6,geo5)
geo8 = setdiff(geo7,geo4)

#using AbstractTrees
#print_tree(stdout,get_tree(geo8))

cutgeo = cut(model,geo8)

trian = Triangulation(model)
writevtk(trian,joinpath(d,"trian"))

trian1 = Triangulation(cutgeo)
writevtk(trian1,joinpath(d,"trian1"))

trian4_Γ = EmbeddedBoundary(cutgeo,geo8,geo4)
writevtk(trian4_Γ,joinpath(d,"trian4_G"))

trian5_Γ = EmbeddedBoundary(cutgeo,geo8,geo5)
writevtk(trian5_Γ,joinpath(d,"trian5_G"))

trian6_Γ = EmbeddedBoundary(cutgeo,geo8,geo6)
writevtk(trian6_Γ,joinpath(d,"trian6_G"))

n = 5
partition = (n,n)
domain = (-0.01,2.01,-0.01,1.01)
model = CartesianDiscreteModel(domain,partition)
trian = Triangulation(model)
#writevtk(trian,"trian")

geo = quadrilateral(x0=Point(0,0),d1=VectorValue(1,0),d2=VectorValue(0,1))
test_geometry(geo)
cutgeo = cut(model,geo)

trian_Ω = Triangulation(cutgeo)
quad_Ω = CellQuadrature(trian_Ω,2)
trian_Γ = EmbeddedBoundary(cutgeo)
n_Γ = get_normal_vector(trian_Γ)
#writevtk(trian_Ω,"trian_O")
#writevtk(trian_Γ,"trian_G",cellfields=["n_g"=>n_Γ])

area_1 = sum(integrate(1,quad_Ω))
tol = 1.0e-9
Area_1 = 1.0
@test Area_1 - area_1 < tol


geo1 = quadrilateral(x0=Point(0,0),d1=VectorValue(1,0),d2=VectorValue(1,1))
test_geometry(geo1)
cutgeo = cut(model,geo1)

trian_Ω = Triangulation(cutgeo)
quad_Ω = CellQuadrature(trian_Ω,2)
trian_Γ = EmbeddedBoundary(cutgeo)
writevtk(trian_Γ,joinpath(d,"trian_G"))

n_Γ = get_normal_vector(trian_Γ)
#writevtk(trian_Ω,"trian_O")
#writevtk(trian_Γ,"trian_G",cellfields=["n_g"=>n_Γ])

area_2 = sum(integrate(1,quad_Ω))
Area_2 = 1.0
@test Area_2 - area_2 < tol

domain = (-0.01,2.01,-0.01,3.01)
model = CartesianDiscreteModel(domain,partition)
trian = Triangulation(model)
#writevtk(trian,"trian")

geo2 = quadrilateral(x0=Point(0,0),d1=VectorValue(1,1),d2=VectorValue(1,2))
test_geometry(geo2)
cutgeo = cut(model,geo2)

trian_Ω = Triangulation(cutgeo)
quad_Ω = CellQuadrature(trian_Ω,2)
trian_Γ = EmbeddedBoundary(cutgeo)
n_Γ = get_normal_vector(trian_Γ)
#writevtk(trian_Ω,"trian_O")
#writevtk(trian_Γ,"trian_G",cellfields=["n_g"=>n_Γ])

area_3 = sum(integrate(1,quad_Ω))
Area_3 = 1.0
@test Area_1 - area_1 < tol
  
end # tmpdir

end # module
