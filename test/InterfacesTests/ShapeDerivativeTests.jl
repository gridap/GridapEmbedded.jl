module ShapeDerivativeTests

using ForwardDiff
using GridapEmbedded
using GridapEmbedded.LevelSetCutters
using Gridap
using Gridap.ReferenceFEs
using Gridap.Arrays

L=1
n = 5
dimensions = 2; domain = (0,L,0,L); cells=(n,n)
bgmodel = CartesianDiscreteModel(domain,cells)
point_to_coords = collect1d(get_node_coordinates(bgmodel))

geoc = disk(L/4,x0=Point(L/2,L/2)) 
geoc = discretize(geoc,bgmodel)
ϕ = geoc.tree.data[1] # signed distance function

geo = DiscreteGeometry(ϕ,point_to_coords,name="")
cutgeo = cut(bgmodel,geo)

Ω = Triangulation(cutgeo)

order = 1
model = DiscreteModel(cutgeo)
V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order),conformity=:H1)
u(x) = x[1]^2 - x[2]
uₕ = interpolate(u,V) 

function f(ϕ)
  ϕ = collect1d(ϕ)
  geo = DiscreteGeometry(ϕ,point_to_coords,name="")
  cutgeo = cut(bgmodel,geo)
  Ω = Triangulation(cutgeo)
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dc = (∫(uₕ)dΩ)
  sum(map(sum,values(dc.dict)))
end

df(ϕ) = ForwardDiff.gradient(f,ϕ) # gradient of an integral wrt the signed distance function

df(ϕ)

#using FiniteDifferences
#using Test
#dfFD(ϕ) = grad(central_fdm(2, 1), f, ϕ)[1]
#tol=1e-8
#@test sum((df(ϕ)) - (dfFD(ϕ))) < tol

end 
