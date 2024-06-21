module PeriodicPoissonAgFEMTests

using Gridap
using GridapEmbedded
using Test

function update_labels!(e::Integer,model::CartesianDiscreteModel,f_Γ::Function,name::String)
  mask = mark_nodes(f_Γ,model)
  _update_labels_locally!(e,model,mask,name)
  nothing
end

function _update_labels_locally!(e,model::CartesianDiscreteModel{2},mask,name)
  topo   = Gridap.Geometry.get_grid_topology(model)
  labels = Gridap.Geometry.get_face_labeling(model)
  cell_to_entity = labels.d_to_dface_to_entity[end]
  entity = maximum(cell_to_entity) + e
  # Vertices
  vtxs_Γ = findall(mask)
  vtx_edge_connectivity = Array(Gridap.Geometry.get_faces(topo,0,1)[vtxs_Γ])
  # Edges
  edge_entries = [findall(x->any(x .∈  vtx_edge_connectivity[1:end.!=j]),
    vtx_edge_connectivity[j]) for j = 1:length(vtx_edge_connectivity)]
  edge_Γ = unique(reduce(vcat,getindex.(vtx_edge_connectivity,edge_entries),init=[]))
  labels.d_to_dface_to_entity[1][vtxs_Γ] .= entity
  labels.d_to_dface_to_entity[2][edge_Γ] .= entity
  add_tag!(labels,name,[entity])
  return cell_to_entity
end

function mark_nodes(f,model::DiscreteModel)
  topo   = Gridap.Geometry.get_grid_topology(model)
  coords = Gridap.Geometry.get_vertex_coordinates(topo)
  mask = map(f,coords)
  return mask
end

u(x) = (x[1]-0.5)^2 + 2(x[2]-0.5)^2
f(x) = -Δ(u)(x)
ud(x) = u(x)

# R = 0.3
# geo1 = square(;L=2)
# geo2 = disk(R,x0=Point(0.5,0.5))

geom = disk(0.55,x0=Point(0.5,0.5))
# geom = setdiff(geo1,geo2)

n = 31
partition = (n,n)
domain = (0,1,0,1)
bgmodel = CartesianDiscreteModel(domain,partition;isperiodic=(true,true))
update_labels!(1,bgmodel,x->iszero(x[1]) || iszero(x[2]),"outer_boundary")

dp = 1
const h = dp/n

cutgeo = cut(bgmodel,geom)

strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutgeo)

Ω_bg = Triangulation(bgmodel)
Ω_act = Triangulation(cutgeo,ACTIVE)
Ω = Triangulation(cutgeo,PHYSICAL)
Γ = EmbeddedBoundary(cutgeo)

n_Γ = get_normal_vector(Γ)

order = 2
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)

model = get_active_model(Ω_act)
Vstd = FESpace(Ω_act,FiniteElements(PhysicalDomain(),model,lagrangian,Float64,order);dirichlet_tags="outer_boundary")

V = AgFEMSpace(Vstd,aggregates)
U = TrialFESpace(V,ud)

const γd = 10.0

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

#colors = color_aggregates(aggregates,bgmodel)
#writevtk(Ω_bg,"trian",celldata=["aggregate"=>aggregates,"color"=>colors],cellfields=["uh"=>uh])
# writevtk(Ω,"trian_O",cellfields=["uh"=>uh,"u"=>u])
# writevtk(Γ,"trian_G")
@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

end # module
