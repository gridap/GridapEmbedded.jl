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
reffe = ReferenceFE(:Lagrangian,Float64,order)
V = FESpace(model,reffe,conformity=:H1)

Random.seed!(1234)
v = FEFunction(V,rand(num_free_dofs(V)))
u = interpolate(x->x[1]+x[2],V)

Ωbg = Triangulation(bgmodel)
Ω = Triangulation(cutgeo)
Γu = EmbeddedBoundary(cutgeo)
Γf = BoundaryTriangulation(cutgeo_facets)
Γ = lazy_append(Γu,Γf)
Λ = SkeletonTriangulation(cutgeo_facets)

dΩ = LebesgueMeasure(Ω,2*order)
dΓ = LebesgueMeasure(Γ,2*order)
dΛ = LebesgueMeasure(Λ,2*order)

n_Γ = get_normal_vector(Γ)
n_Λ = get_normal_vector(Λ)

# Check divergence theorem
a = sum( ∫( ∇(v)⋅∇(u) )*dΩ )
b = sum( ∫( v*n_Γ⋅∇(u) )*dΓ )
@test abs(a-b) < 1.0e-9

a = sum( ∫( jump(u) )*dΛ )
@test abs(a) < 1.0e-9

a = sum( ∫( jump(v) )*dΛ )
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
#writevtk(Ωbg,"trian")
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
