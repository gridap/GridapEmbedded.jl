module PeriodicDistributedDiscreteGeometryPoissonTest

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

using GridapDistributed: DistributedDiscreteModel

using Gridap.Geometry
using Gridap.Geometry: get_vertex_coordinates,get_faces

using GridapEmbedded.CSG
using GridapEmbedded.LevelSetCutters

function update_labels!(e::Integer,model::DistributedDiscreteModel,f_Γ::Function,name::String)
  mask = mark_nodes(f_Γ,model)
  cell_to_entity = map(local_views(model),local_views(mask)) do model,mask
    _update_labels_locally!(e,model,mask,name)
  end
  cell_gids=get_cell_gids(model)
  cache=GridapDistributed.fetch_vector_ghost_values_cache(cell_to_entity,partition(cell_gids))
  GridapDistributed.fetch_vector_ghost_values!(cell_to_entity,cache)
  nothing
end

function _update_labels_locally!(e,model::CartesianDiscreteModel{2},mask,name)
  topo   = get_grid_topology(model)
  labels = get_face_labeling(model)
  cell_to_entity = labels.d_to_dface_to_entity[end]
  entity = maximum(cell_to_entity) + e
  # Vertices
  vtxs_Γ = findall(mask)
  vtx_edge_connectivity = Array(get_faces(topo,0,1)[vtxs_Γ])
  # Edges
  edge_entries = [findall(x->any(x .∈  vtx_edge_connectivity[1:end.!=j]),
    vtx_edge_connectivity[j]) for j = 1:length(vtx_edge_connectivity)]
  edge_Γ = unique(reduce(vcat,getindex.(vtx_edge_connectivity,edge_entries),init=[]))
  labels.d_to_dface_to_entity[1][vtxs_Γ] .= entity
  labels.d_to_dface_to_entity[2][edge_Γ] .= entity
  add_tag!(labels,name,[entity])
  return cell_to_entity
end

function mark_nodes(f,model::DistributedDiscreteModel)
  local_masks = map(local_views(model)) do model
    mark_nodes(f,model)
  end
  gids = get_face_gids(model,0)
  mask = PVector(local_masks,partition(gids))
  assemble!(|,mask) |> fetch  # Ghosts -> Owned with `or` applied
  consistent!(mask) |> fetch  # Owned  -> Ghost
  return mask
end

function mark_nodes(f,model::DiscreteModel)
  topo   = get_grid_topology(model)
  coords = get_vertex_coordinates(topo)
  mask = map(f,coords)
  return mask
end

function main(distribute,parts;
  threshold=1,
  n=8,
  cells=(n,n),
  geometry=:circle)

  ranks = distribute(LinearIndices((prod(parts),)))

  u(x) = (x[1]-0.5)^2 + (x[2]-0.5)^2
  f(x) = -Δ(u)(x)
  ud(x) = u(x)

  geometries = Dict(
    :circle => circle_geometry,
    :remotes => remotes_geometry,
  )

  bgmodel,_geo = geometries[geometry](ranks,parts,cells)
  update_labels!(1,bgmodel,x->iszero(x[1]) || iszero(x[2]),"outer_boundary")

  geo = discretize(_geo,bgmodel)

  D = 2
  cell_meas = map(get_cell_measure∘Triangulation,local_views(bgmodel))
  meas = map(first,cell_meas) |> PartitionedArrays.getany
  h = meas^(1/D)

  cutgeo = cut(bgmodel,geo)
  cutgeo_facets = cut_facets(bgmodel,geo)

  strategy = AggregateCutCellsByThreshold(threshold)
  bgmodel,cutgeo,aggregates = aggregate(strategy,cutgeo)

  Ω_bg = Triangulation(bgmodel)
  Ω_act = Triangulation(cutgeo,ACTIVE)
  Ω = Triangulation(cutgeo,PHYSICAL)
  Γ = EmbeddedBoundary(cutgeo)

  n_Γ = get_normal_vector(Γ)

  order = 2
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓ = Measure(Γ,degree)

  reffe = ReferenceFE(lagrangian,Float64,order)

  Vstd = FESpace(Ω_act,reffe;dirichlet_tags="outer_boundary")

  V = AgFEMSpace(bgmodel,Vstd,aggregates)
  U = TrialFESpace(V,ud)

  γd = 10.0

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

  #
  colors = map(color_aggregates,aggregates,local_views(bgmodel))
  gids = get_cell_gids(bgmodel)


  global_aggregates = map(aggregates,local_to_global(gids)) do agg,gid
    map(i-> i==0 ? 0 : gid[i],agg)
  end
  own_aggregates = map(global_aggregates,own_to_local(gids)) do agg,oid
    map(Reindex(agg),oid)
  end
  own_colors = map(colors,own_to_local(gids)) do col,oid
    map(Reindex(col),oid)
  end

  path = mktempdir()
  writevtk(Ω_bg,joinpath(path,"trian"),
    celldata=[
      "aggregate"=>own_aggregates,
      "color"=>own_colors,
      "gid"=>own_to_global(gids)])#,
  #  cellfields=["uh"=>uh])

  writevtk(Ω,joinpath(path,"trian_O"),cellfields=["uh"=>uh])
  writevtk(Γ,joinpath(path,"trian_G"))
  @test el2/ul2 < 1.e-8
  @test eh1/uh1 < 1.e-7

end

function circle_geometry(ranks,parts,cells)
  p0 = Point(0.5,0.5)
  R = 0.55
  geo = disk(R,x0=p0)
  bgmodel = CartesianDiscreteModel(ranks,parts,(0,1,0,1),cells;isperiodic=(true,true))
  bgmodel,geo
end

function remotes_geometry(ranks,parts,cells)
  x0 = Point(0.0,0.4)
  d1 = VectorValue(1.0,0.0)
  d2 = VectorValue(0.0,0.2)
  geo1 = quadrilateral(;x0=x0,d1=d1,d2=d2)

  x0 = Point(0.4,0.0)
  d1 = VectorValue(0.2,0.0)
  d2 = VectorValue(0.0,1.0)
  geo2 = quadrilateral(;x0=x0,d1=d1,d2=d2)
  geo = union(geo1,geo2)

  domain = (0, 1, 0, 1)
  bgmodel = CartesianDiscreteModel(ranks,parts,domain,cells;isperiodic=(true,true))
  bgmodel,geo
end

end # module