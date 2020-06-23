module EmbeddedFacetDiscretizationsTests

using Random
using Test
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters

n = 10
partition = (n,n)
domain = (0,1,0,1)

R = 0.72
geo = disk(R,x0=Point(1,1))

bgmodel = CartesianDiscreteModel(domain,partition)

cutgeo = cut(bgmodel,geo)
cutgeo_facets = cut_facets(bgmodel,geo)

model = DiscreteModel(cutgeo)

order = 1
V = FESpace(
  model=model,conformity=:H1,reffe=:Lagrangian,order=order,valuetype=Float64)

U = TrialFESpace(V)

Random.seed!(1234)
v = FEFunction(V,rand(num_free_dofs(V)))

u = interpolate(U,x->x[1]+x[2])

trian = Triangulation(bgmodel)
trian_Ω = Triangulation(cutgeo)
trian_Γu = EmbeddedBoundary(cutgeo)
trian_Γf = BoundaryTriangulation(cutgeo_facets)
trian_Γ = lazy_append(trian_Γu,trian_Γf)
trian_sΩ = SkeletonTriangulation(cutgeo_facets)

quad_Ω = CellQuadrature(trian_Ω,2*order)
quad_Γ = CellQuadrature(trian_Γ,2*order)
quad_sΩ = CellQuadrature(trian_sΩ,2*order)

n_Γ = get_normal_vector(trian_Γ)
n_sΩ = get_normal_vector(trian_sΩ)

v_Ω = restrict(v,trian_Ω)
v_Γ = restrict(v,trian_Γ)
u_Ω = restrict(u,trian_Ω)
u_Γ = restrict(u,trian_Γ)

v_sΩ = restrict(v,trian_sΩ)
u_sΩ = restrict(u,trian_sΩ)

# Check divergence theorem
a = sum( integrate(∇(v_Ω)⋅∇(u_Ω),trian_Ω,quad_Ω) )
b = sum( integrate(v_Γ*n_Γ⋅∇(u_Γ),trian_Γ,quad_Γ) )
@test abs(a-b) < 1.0e-9

a = sum(integrate(jump(u_sΩ),trian_sΩ,quad_sΩ))
@test abs(a) < 1.0e-9

a = sum(integrate(jump(v_sΩ),trian_sΩ,quad_sΩ))
@test abs(a) < 1.0e-9

#celldata_Ω = ["bgcell"=>collect(Int,get_cell_id(trian_Ω))]
#celldata_Γ = ["bgcell"=>collect(Int,get_cell_id(trian_Γ))]
#cellfields_Ω = ["v"=>v_Ω,"u"=>u_Ω]
#cellfields_Γ = ["normal"=>n_Γ,"v"=>v_Γ,"u"=>u_Γ]
#
#celldata_sΩ = [
#  "bgcell_left"=>collect(Int,get_cell_id(trian_sΩ).left),
#  "bgcell_right"=>collect(Int,get_cell_id(trian_sΩ).right)]
#cellfields_sΩ = ["normal"=> n_sΩ,"jump_v"=>jump(v_sΩ),"jump_u"=>jump(u_sΩ)]
#
#writevtk(trian,"trian")
#writevtk(trian_Ω,"trian_O",celldata=celldata_Ω,cellfields=cellfields_Ω)
#writevtk(trian_Γ,"trian_G",celldata=celldata_Γ,cellfields=cellfields_Γ)
#writevtk(trian_sΩ,"trian_sO",celldata=celldata_sΩ,cellfields=cellfields_sΩ)

n = 10
partition = (n,n,n)
domain = (0,1,0,1,0,1)

R = 0.49
geo = disk(R,x0=Point(.5,.5,.5))

bgmodel = CartesianDiscreteModel(domain,partition)

cutgeo = cut(bgmodel,geo)
cutgeo_facets = cut_facets(bgmodel,geo)

trian_s = SkeletonTriangulation(bgmodel)
trian_sΩ = SkeletonTriangulation(cutgeo_facets,trian_s,geo,(CUTIN,IN))
trian_sΩo = SkeletonTriangulation(cutgeo_facets,trian_s,geo,(CUTOUT,OUT))

#writevtk(trian_s,"trian_s")
#writevtk(trian_sΩ,"trian_sO")
#writevtk(trian_sΩo,"trian_sOo")

end # module
