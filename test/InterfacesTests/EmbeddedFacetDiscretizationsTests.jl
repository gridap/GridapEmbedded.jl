module EmbeddedFacetDiscretizationsTests

using Random
using Test
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters

##########################################
# 2D tests 

n = 10
partition = (n,n)
domain = (0,1,0,1)

R = 0.72
geo = disk(R,x0=Point(1,1))

bgmodel = CartesianDiscreteModel(domain,partition)

cutgeo = cut(bgmodel,geo)
cutgeo_facets = cut_facets(bgmodel,geo)

Ωact = Triangulation(cutgeo,ACTIVE)

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(Ωact,reffe)

Random.seed!(1234)
v = FEFunction(V,rand(num_free_dofs(V)))
u = interpolate(x->x[1]+x[2],V)

Ωbg = Triangulation(bgmodel)
Ω = Triangulation(cutgeo,PHYSICAL)
Γu = EmbeddedBoundary(cutgeo)
Γf = BoundaryTriangulation(cutgeo_facets,PHYSICAL)
Γ = lazy_append(Γu,Γf)
Λ = SkeletonTriangulation(cutgeo_facets,PHYSICAL)

face_model = get_active_model(Γu)
Σb = BoundaryTriangulation(Γu)
Σi = SkeletonTriangulation(Γu)

test_triangulation(Ω)
test_triangulation(Γ)
test_triangulation(Λ)

dΩ = Measure(Ω,2*order)
dΓ = Measure(Γ,2*order)
dΛ = Measure(Λ,2*order)

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

D = num_cell_dims(bgmodel)
celldata_Ω = ["bgcell"=>collect(Int,get_glue(Ω,Val(D)).tface_to_mface)]
celldata_Γ = ["bgcell"=>collect(Int,get_glue(Γ,Val(D)).tface_to_mface)]
cellfields_Ω = ["v"=>v,"u"=>u]
cellfields_Γ = ["normal"=>n_Γ,"v"=>v,"u"=>u]
celldata_Λ = [
 "bgcell_left"=>collect(Int,get_glue(Λ.⁺,Val(D)).tface_to_mface),
 "bgcell_right"=>collect(Int,get_glue(Λ.⁻,Val(D)).tface_to_mface)]
cellfields_Λ = ["normal"=> n_Λ.⁺,"jump_v"=>jump(v),"jump_u"=>jump(u)]

d = mktempdir()
try
  writevtk(Ωbg,joinpath(d,"trian"),append=false)
  writevtk(Ω,joinpath(d,"trian_O"),celldata=celldata_Ω,cellfields=cellfields_Ω,append=false)
  writevtk(Γ,joinpath(d,"trian_G"),celldata=celldata_Γ,cellfields=cellfields_Γ,append=false)
  writevtk(Λ,joinpath(d,"trian_sO"),celldata=celldata_Λ,cellfields=cellfields_Λ,append=false)
  writevtk(Γf,joinpath(d,"trian_Gf"),append=false)
finally
  rm(d,recursive=true)
end

##########################################
# 3D tests 

n = 10
partition = (n,n,n)
domain = (0,1,0,1,0,1)

R = 0.49
geo = disk(R,x0=Point(.5,.5,.5))

bgmodel = CartesianDiscreteModel(domain,partition)

cutgeo = cut(bgmodel,geo)
cutgeo_facets = cut_facets(bgmodel,geo)

trian_s = SkeletonTriangulation(bgmodel)
trian_sΩ = SkeletonTriangulation(trian_s,cutgeo_facets,PHYSICAL_IN,geo)
trian_sΩo = SkeletonTriangulation(trian_s,cutgeo_facets,PHYSICAL_OUT,geo)

Γu = EmbeddedBoundary(cutgeo)
face_model = get_active_model(Γu)
Σb = BoundaryTriangulation(Γu)
Σi = SkeletonTriangulation(Γu)

d = mktempdir()
try
  writevtk(trian_s,joinpath(d,"trian_s"))
  writevtk(trian_sΩ,joinpath(d,"trian_sO"))
  writevtk(trian_sΩo,joinpath(d,"trian_sOo"))
finally
  rm(d,recursive=true)
end

end # module
