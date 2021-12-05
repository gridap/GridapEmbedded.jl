module ShapeDiff

#using ForwardDiff
using ReverseDiff
using GridapEmbedded
using GridapEmbedded.LevelSetCutters
using Gridap
using Gridap.ReferenceFEs
using Gridap.Arrays

domain = (0,1,0,1); cells=(10,10)
bgmodel = CartesianDiscreteModel(domain,cells)

point_to_coords = collect1d(get_node_coordinates(bgmodel))
geo1 = disk(0.3,x0=Point(0.5,0.5))
geo1_d = discretize(geo1,bgmodel) 
lvl_set = geo1_d.tree.data[1]

geo = DiscreteGeometry(lvl_set,point_to_coords,name="")
order=1

function residual(ϕₘ)
  lvl_set = collect1d(ϕₘ)
  geo = DiscreteGeometry(lvl_set,point_to_coords,name="")
  cutgeo = cut(bgmodel,geo)
  Ω = Triangulation(cutgeo,PHYSICAL)
  Ω1_act = Triangulation(cutgeo,ACTIVE)
  dΩ = Measure(Ω,2)
  Vstd = FESpace(Ω1_act,ReferenceFE(lagrangian,Float64,order),conformity=:H1)  
  u0=interpolate(1,Vstd)
  sum(∫(u0)dΩ)
end

residual(lvl_set)

ϕ = lvl_set
ϕₘ =  reshape(ϕ,(length(ϕ),1))
#dϕₘ = ForwardDiff.gradient(residual,ϕₘ)
dϕₘ = ReverseDiff.gradient(residual,ϕₘ)

end 