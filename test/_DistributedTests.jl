module DistributedPoissonTests

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

using GridapDistributed: DistributedDiscreteModel
using GridapDistributed: DistributedGridapType
using GridapDistributed: DistributedTriangulation
using GridapDistributed: DistributedFESpace
using GridapDistributed: generate_gids

import Gridap.Geometry: Triangulation
import Gridap.Geometry: get_background_model
import GridapEmbedded.Interfaces: cut
import GridapEmbedded.Interfaces: Cutter
import GridapEmbedded.Interfaces: EmbeddedBoundary
import GridapEmbedded.AgFEM: aggregate
import GridapEmbedded.AgFEM: AgFEMSpace
import GridapDistributed: local_views
import GridapDistributed: remove_ghost_cells
using GridapEmbedded.Interfaces: SubFacetTriangulation
using Gridap.Geometry: AppendedTriangulation
using Gridap.Arrays

struct DistributedEmbeddedDiscretization{Dp,T,A,B} <: DistributedGridapType
  discretizations::A
  model::B
  function DistributedEmbeddedDiscretization(
    d::AbstractPData{<:EmbeddedDiscretization{Dp,T}},
    model::DistributedDiscreteModel) where {Dp,T}
    A = typeof(d)
    B = typeof(model)
    new{Dp,T,A,B}(d,model)
  end
end

local_views(a::DistributedEmbeddedDiscretization) = a.discretizations

get_background_model(a::DistributedEmbeddedDiscretization) = a.model

function cut(cutter::Cutter,bgmodel::DistributedDiscreteModel,args...)
  discretizations = map_parts(local_views(bgmodel)) do lmodel
    cut(cutter,lmodel,args...)
  end
  DistributedEmbeddedDiscretization(discretizations,bgmodel)
end

function cut(bgmodel::DistributedDiscreteModel,args...)
  cutter = LevelSetCutter()
  cut(cutter,bgmodel,args...)
end

function Triangulation(cutgeo::DistributedEmbeddedDiscretization,args...)
  trians = map_parts(local_views(cutgeo)) do lcutgeo
    Triangulation(lcutgeo,args...)
  end
  bgmodel = get_background_model(cutgeo)
  DistributedTriangulation(trians,bgmodel)
end

function EmbeddedBoundary(cutgeo::DistributedEmbeddedDiscretization,args...)
  trians = map_parts(local_views(cutgeo)) do lcutgeo
    EmbeddedBoundary(lcutgeo,args...)
  end
  bgmodel = get_background_model(cutgeo)
  DistributedTriangulation(trians,bgmodel)
end

function aggregate(stragegy,cutgeo::DistributedEmbeddedDiscretization)
  map_parts(local_views(cutgeo)) do lcutgeo
    aggregate(stragegy,lcutgeo)
  end
  # TODO: parallel aggregation here
end

function AgFEMSpace(
    bgmodel::DistributedDiscreteModel,
    f::DistributedFESpace,
    bgcell_to_bgcellin::SequentialData{<:AbstractVector},
    g::DistributedFESpace=f)

  fagg = map_parts(AgFEMSpace,local_views(f),bgcell_to_bgcellin,local_views(g))
  gids = generate_gids(bgmodel,fagg)
  vector_type = get_vector_type(f)
  DistributedSingleFieldFESpace(fagg,gids,vector_type)
end

function global_aggregates(bgmodel::DistributedDiscreteModel,aggregates)
  gids = map_parts(get_lid_to_gid,local_views(get_cell_gids(bgmodel)))
  map_parts(gids,aggregates) do gid,agg
    map( i ->i == 0 ? 0 : gid[i], agg)
  end
end

function remove_ghost_cells(trian::DistributedTriangulation)
  model = get_background_model(trian)
  gids = get_cell_gids(model)
  trians = map_parts(local_views(trian),local_views(gids)) do trian,gids
    remove_ghost_cells(trian,gids)
  end
  DistributedTriangulation(trians,model)
end

function remove_ghost_cells(trian::AppendedTriangulation,gids)
  a = remove_ghost_cells(trian.a,gids)
  b = remove_ghost_cells(trian.b,gids)
  lazy_append(a,b)
end

function remove_ghost_cells(trian::SubFacetTriangulation,gids)
  model = get_background_model(trian)
  D     = num_cell_dims(model)
  glue  = get_glue(trian,Val{D}())
  remove_ghost_cells(glue,trian,gids)
end

partition = (2,2)

parts = get_part_ids(SequentialBackend(),partition)

u(x) = x[1] - x[2]
f(x) = -Δ(u)(x)
ud(x) = u(x)

R = 0.5
L = 0.8*(2*R)
p0 = Point(0.0,0.0)

geo = disk(R,x0=p0)

t = 1.1
pmin = p0-t*R
pmax = p0+t*R

n = 10
mesh_partition = (n,n)
bgmodel = CartesianDiscreteModel(parts,pmin,pmax,mesh_partition)
bgtrian = Triangulation(bgmodel)
writevtk(bgtrian,"bgtrian")

dp = pmax - pmin
h = dp[1]/n

cutgeo = cut(bgmodel,geo)

strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutgeo)

Ω_bg = Triangulation(bgmodel)
Ω_act = Triangulation(cutgeo,ACTIVE)
Ω = Triangulation(cutgeo,PHYSICAL)
Γ = EmbeddedBoundary(cutgeo)

Γ = remove_ghost_cells(Γ)
Ω = remove_ghost_cells(Ω)

t1 = get_part(local_views(Ω),1)
writevtk(t1,"t1")

writevtk(Ω,"trian")
writevtk(Ω_act,"trian_act")
n_Γ = get_normal_vector(Γ)

order = 1
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)

reffe = ReferenceFE(lagrangian,Float64,order)

Vstd = FESpace(Ω_act,reffe)

#V = AgFEMSpace(bgmodel,Vstd,aggregates)
V = Vstd
U = TrialFESpace(V)


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

#
colors = map_parts(color_aggregates,aggregates,local_views(bgmodel))
cell_gids = get_cell_gids(bgmodel)
ohids = map_parts(get_lid_to_ohid,local_views(cell_gids))
gids = map_parts(get_lid_to_gid,local_views(cell_gids))
oids = map_parts(i->findall(>(0),i),ohids)

_aggregates = map_parts(aggregates,gids) do agg,gid
  map(i-> i==0 ? 0 : gid[i],agg)
end
_aggregates = map_parts(_aggregates,oids) do agg,oid
  map(Reindex(agg),oid)
end
_colors = map_parts(colors,oids) do col,oid
  map(Reindex(col),oid)
end
_gids = map_parts(gids,oids) do gid,oid
  map(Reindex(gid),oid)
end

writevtk(Ω_bg,"trian",
  celldata=[
    "aggregate"=>_aggregates,
    "color"=>_colors,
    "gid"=>_gids],
  cellfields=["uh"=>uh])

writevtk(Ω,"trian_O",cellfields=["uh"=>uh])
writevtk(Γ,"trian_G")
@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7

end # module
