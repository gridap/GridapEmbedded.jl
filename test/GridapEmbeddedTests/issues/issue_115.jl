using Gridap, GridapEmbedded
using Gridap.Geometry, Gridap.Adaptivity
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.LevelSetCutters: SubCellTriangulation
using Test

### Test 1
n = 12
base_model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(n,n)))
ref_model = refine(base_model, refinement_method = "barycentric")
model = get_model(ref_model)
φ1((x,y)) = (x-0.5)^2 + (y-0.5)^2 - 0.3^2
φ2((x,y)) = (x-0.5)^2 + (y-0.5)^2 - (0.3-0.5/n)^2
V_φ = TestFESpace(model,ReferenceFE(lagrangian,Float64,1))
φh1 = interpolate(φ1,V_φ); φh2 = interpolate(φ2,V_φ)

geo1 = DiscreteGeometry(φh1,model,name="φ1")
geo2 = DiscreteGeometry(φh2,model,name="φ2")
Ω1_geo = setdiff(!geo1,geo2,name="Ω1")
Ω2_geo = setdiff(geo2,geo1,name="Ω2")
cutgeo = cut(model,union(Ω1_geo,Ω2_geo))

Ω2 = Triangulation(cutgeo,PHYSICAL,"Ω2")
Ω2_act = Triangulation(cutgeo,ACTIVE,"Ω2")
Λ_Ω2 = GhostSkeleton(cutgeo,"Ω2")

@test iszero(num_cells(Ω2))
@test iszero(num_cells(Ω2_act))
@test iszero(num_cells(Λ_Ω2))

### Test 2
n = 48
base_model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(n,n)))
ref_model = refine(base_model, refinement_method = "barycentric")
model = get_model(ref_model)
f((x,y),a1,a2,b1,b2,c) = cos(a1*pi*(x-a2))*cos(b1*pi*(y-b2)) - c
φ1((x,y)) = f((x,y),4,1/8,4,1/8,0.5)
φ2((x,y)) = f((x,y),4,3/8,4,1/8,0.5)
φ3((x,y)) = f((x,y),4,3/8,4,1/8,0.5-0.5/n)
V_φ = TestFESpace(model,ReferenceFE(lagrangian,Float64,1))
φh1 = interpolate(φ1,V_φ); φh2 = interpolate(φ2,V_φ); φh3 = interpolate(φ3,V_φ)

geo1 = DiscreteGeometry(φh1,model,name="φ1")
geo2 = DiscreteGeometry(φh2,model,name="φ2")
geo3 = DiscreteGeometry(φh3,model,name="φ3")
Ω1_geo = intersect(geo1,geo2,name="Ω1")
Ω2_geo = setdiff(geo3,Ω1_geo,name="Ω2")
cutgeo = cut(model,union(Ω1_geo,Ω2_geo))

Ω2 = Triangulation(cutgeo,PHYSICAL,"Ω2")
Ω2_act = Triangulation(cutgeo,ACTIVE,"Ω2")
Λ_Ω2 = GhostSkeleton(cutgeo,"Ω2")

# Check cells are in subcell to bgcell list or in cell to bgcell list
for c in Ω2_act.tface_to_mface
  @test c ∈ Ω2.a.subcells.cell_to_bgcell || c ∈ Ω2.b.tface_to_mface
end
for c in Λ_Ω2.plus.glue.face_to_cell
  @test c ∈ Ω2.a.subcells.cell_to_bgcell || c ∈ Ω2.b.tface_to_mface
end
for c in Λ_Ω2.minus.glue.face_to_cell
  @test c ∈ Ω2.a.subcells.cell_to_bgcell || c ∈ Ω2.b.tface_to_mface
end

writevtk(Triangulation(model),"Omega")
writevtk(Triangulation(cutgeo,"Ω2"),"Omega2")
writevtk(Triangulation(cutgeo,ACTIVE,"Ω2"),"Omega2_act")
writevtk(Λ_Ω2,"Omega2_ghost_skel")