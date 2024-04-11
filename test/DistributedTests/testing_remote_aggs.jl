module testing_remote_aggs

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

distribute = PartitionedArrays.DebugArray


np = (3,1)
ranks = distribute(LinearIndices((prod(np),)))

nc = (9,9)

L = 1

# C shape
x0 = Point(0.1,0.1)
d1 = VectorValue(0.8,0.0)
d2 = VectorValue(0.0,0.05)
geo1 = quadrilateral(;x0=x0,d1=d1,d2=d2)

x0 = Point(0.15,0.125)
d1 = VectorValue(0.2,0.0)
d2 = VectorValue(0.0,0.6)
geo2 = quadrilateral(;x0=x0,d1=d1,d2=d2)

x0 = Point(0.1,0.6)
d1 = VectorValue(0.8,0.0)
d2 = VectorValue(0.0,0.3)
geo3 = quadrilateral(;x0=x0,d1=d1,d2=d2)

geo12 = union(geo1,geo2)
geo = union(geo12,geo3)
geo = geo12

# Cirlce (test)
R = 0.35
p0 = Point(0.5,0.5)
geo = disk(R,x0=p0)
geo = geo3


# Bigger  cut cells (test)
x0 = Point(0.05,0.05)
d1 = VectorValue(0.9,0.0)
d2 = VectorValue(0.0,0.15)
geo1 = quadrilateral(;x0=x0,d1=d1,d2=d2)

x0 = Point(0.15,0.125)
d1 = VectorValue(0.25,0.0)
d2 = VectorValue(0.0,0.6)
geo2 = quadrilateral(;x0=x0,d1=d1,d2=d2)

x0 = Point(0.05,0.6)
d1 = VectorValue(0.9,0.0)
d2 = VectorValue(0.0,0.35)
geo3 = quadrilateral(;x0=x0,d1=d1,d2=d2)

geo12 = union(geo1,geo2)
geo = union(geo12,geo3)

geo = geo12

domain = (0, 1, 0, 1)
bgmodel = CartesianDiscreteModel(ranks,np,domain,nc)

cutgeo = cut(bgmodel,geo)

strategy = AggregateCutCellsByThreshold(1.0)


bgmodel,cutgeo,aggregates = aggregate(strategy,cutgeo,geo,IN);


Ω_bg = Triangulation(bgmodel)
Ω_act = Triangulation(cutgeo,ACTIVE)
Ω = Triangulation(cutgeo,PHYSICAL)
Γ = EmbeddedBoundary(cutgeo)
n_Γ = get_normal_vector(Γ)

order = 1
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)

reffe = ReferenceFE(lagrangian,Float64,order)

Vstd = FESpace(Ω_act,reffe)
V = AgFEMSpace(bgmodel,Vstd,aggregates)

U = TrialFESpace(V)


u(x) = x[1] - x[2]
f(x) = -Δ(u)(x)
ud(x) = u(x)

γd = 10.0
h = 1/nc[1]

a(u,v) =
  ∫( ∇(v)⋅∇(u) ) * dΩ +
  ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u ) * dΓ

l(v) =
  ∫( v*f ) * dΩ +
  ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud ) * dΓ

op = AffineFEOperator(a,l,U,V)
uh = solve(op)


e = u - uh

l2(u) = sqrt(sum( ∫( u*u )*dΩ ))
h1(u) = sqrt(sum( ∫( u*u + ∇(u)⋅∇(u) )*dΩ ))

el2 = l2(e)
eh1 = h1(e)
ul2 = l2(uh)
uh1 = h1(uh)


Ω = Triangulation(cutgeo)
Γ = EmbeddedBoundary(cutgeo)
Ωbg = Triangulation(bgmodel)
Ω_act = Triangulation(cutgeo,ACTIVE)

writevtk(Ω,"trian");
writevtk(Γ,"bnd");
writevtk(Ωbg,"bg_trian");
writevtk(Ω_act,"act_trian");

writevtk(Ω,"trian",
  cellfields=["uh"=>uh,"u"=>u,"e"=>e],);



# # Debug VTK
# map(local_views(_bgmodel),ranks) do m,p
#   writevtk(get_grid(m),"bgmodel_$p")
# end


# map(local_views(uh),local_views(bgmodel),ranks) do uh,m,p
#   trian = Triangulation(m)
#   writevtk(trian,"ltrian_$p",cellfields=["uh"=>uh])
# end


# map(local_views(cutgeo),ranks) do c,p
#   m = get_background_model(c)
#   c_io = compute_bgcell_to_inoutcut(c,c.geo)
#   writevtk(Triangulation(m),"bgtrian_$p",celldata=["inout"=>c_io])
# end


# TODO:
# add geometries
# add tests with mpi
# PR
end # module
