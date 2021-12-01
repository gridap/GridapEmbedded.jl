module AgFEMSpacesTests

using Test
using Gridap
using GridapEmbedded
using Gridap.Geometry: get_active_model

const R = 1.4
geom = disk(R,x0=Point(1.0,1.7))
n = 3
partition = (n,n)

domain = (0,1,0,1)
bgmodel = CartesianDiscreteModel(domain,partition)

cutdisc = cut(bgmodel,geom)

Ω_bg = Triangulation(bgmodel)
Ω_ac = Triangulation(cutdisc,ACTIVE)
Ω = Triangulation(cutdisc,PHYSICAL)
Ω_in = Triangulation(cutdisc,IN)

dΩ_bg = Measure(Ω_bg,2)
dΩ = Measure(Ω,2)
dΩ_in = Measure(Ω_in,2)

model = get_active_model(Ω_ac)
order = 1

# In the physical domain

cell_fe = FiniteElements(PhysicalDomain(),model,lagrangian,Float64,order)
V = FESpace(Ω_ac,cell_fe)

cell_to_cellin = [0,0,9,8,8,9,8,8,9]
Vagg = AgFEMSpace(V,cell_to_cellin)

v(x) = x[1] + 2*x[2]
vhagg = interpolate(v,Vagg)

tol = 10e-9
@test sum( ∫(abs2(v-vhagg))dΩ ) < tol
@test sum( ∫(abs2(v-vhagg))dΩ_in ) < tol

vh = FEFunction(V,rand(num_free_dofs(V)))
vhagg = interpolate(vh,Vagg)
@test sum( ∫(abs2(vh-vhagg))dΩ_in ) < tol

# In the reference space

reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(Ω_ac,reffe)
cell_to_cellin = [0,0,9,8,8,9,8,8,9]
Vagg = AgFEMSpace(V,cell_to_cellin)

v(x) = x[1] + 2*x[2]
vhagg = interpolate(v,Vagg)

tol = 10e-9
@test sum( ∫(abs2(v-vhagg))dΩ ) < tol
@test sum( ∫(abs2(v-vhagg))dΩ_in ) < tol

vh = FEFunction(V,rand(num_free_dofs(V)))
vhagg = interpolate(vh,Vagg)
@test sum( ∫(abs2(vh-vhagg))dΩ_in ) < tol

#cellfields = ["vh"=>vh,"vhagg"=>vhagg,"e"=>vh-vhagg]
#writevtk(Ω_bg,"trian_bg",nsubcells=10,cellfields=cellfields)
#writevtk(Ω_in,"trian_in",nsubcells=10,cellfields=cellfields)
#writevtk(Ω,"trian_phys",cellfields=cellfields)

end # module
