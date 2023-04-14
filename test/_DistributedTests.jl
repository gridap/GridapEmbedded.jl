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
using GridapDistributed: DistributedSingleFieldFESpace
using GridapDistributed: generate_gids
using GridapDistributed: generate_cell_gids
using GridapDistributed: _find_vector_type
using GridapEmbedded.Interfaces: SubFacetTriangulation
using Gridap.Geometry
using Gridap.Geometry: AppendedTriangulation
using Gridap.Geometry: FaceToFaceGlue
using Gridap.FESpaces: SingleFieldFESpace
using Gridap.Arrays
using Gridap.Helpers
using PartitionedArrays: LinearGidToPart
using PartitionedArrays: get_part_ids
using GridapEmbedded.Interfaces: SubCellData
using GridapEmbedded.Interfaces: SubFacetData
using GridapEmbedded.CSG

import Gridap.Geometry: Triangulation
import Gridap.Geometry: get_background_model
import GridapEmbedded.Interfaces: cut
import GridapEmbedded.Interfaces: Cutter
import GridapEmbedded.Interfaces: EmbeddedBoundary
import GridapEmbedded.AgFEM: aggregate
import GridapEmbedded.AgFEM: AgFEMSpace
import GridapDistributed: local_views
import GridapDistributed: remove_ghost_cells
import GridapDistributed: dof_wise_to_cell_wise!
import GridapDistributed: dof_wise_to_cell_wise
import GridapDistributed: generate_gids
import GridapDistributed: add_ghost_cells
import PartitionedArrays: _part_to_firstgid
import GridapEmbedded.Interfaces: compute_bgfacet_to_inoutcut

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
  # TODO: do not compute ghost cells
  #  - cut own_bgmodels
  #  - set (local) bgmodels
  #  - reindex bgcell/facet data
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
  bgcell_to_bgcellin::AbstractPData{<:AbstractVector},
  g::DistributedFESpace=f)

  bgmodel_gids = get_cell_gids(bgmodel)
  spaces = map_parts(
    local_views(f),
    bgcell_to_bgcellin,
    local_views(g),
    local_views(bgmodel_gids)) do f,bgcell_to_bgcellin,g,gids
    AgFEMSpace(f,bgcell_to_bgcellin,g,get_lid_to_gid(gids))
  end
  trians = map_parts(get_triangulation,local_views(f))
  trian = DistributedTriangulation(trians,bgmodel)
  trian = add_ghost_cells(trian)
  trian_gids = generate_cell_gids(trian)
  cell_to_ldofs = map_parts(get_cell_dof_ids,spaces)
  cell_to_ldofs = map_parts(i->map(sort,i),cell_to_ldofs)
  nldofs = map_parts(num_free_dofs,spaces)
  gids = generate_gids(trian_gids,cell_to_ldofs,nldofs)
  gids = add_gid_to_part(gids)
  vector_type = _find_vector_type(spaces,gids)
  DistributedSingleFieldFESpace(spaces,gids,vector_type)
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

## PartitionedArrays.jl

function _part_to_firstgid(ngids::Vector,np::Integer)
  @check length(ngids) == np
  s = zeros(Int32,np)
  s[2:np] = cumsum(ngids)[1:np-1]
  s .+= 1
end

function _gid_to_part(gids::PRange)
  if gids.gid_to_part !== nothing
    gids.gid_to_part
  else
    np = num_parts(gids)
    noids = map_parts(num_oids,local_views(gids))
    ngids = gather_all(noids)
    gid_to_part = map_parts(ngids) do ngids
      part_to_firstgid = _part_to_firstgid(ngids,np)
      LinearGidToPart(sum(ngids),part_to_firstgid)
    end
    map_parts(_test_gid_to_part,local_views(gids),gid_to_part)
    gid_to_part
  end
end

function _test_gid_to_part(i::IndexSet,gid_to_part::AbstractVector)
  for (gid,part) in zip(i.lid_to_gid,i.lid_to_part)
    @check gid_to_part[gid] == part
  end
end

function add_gid_to_part(gids::PRange,gid_to_part=_gid_to_part(gids))
  PRange(
    gids.ngids,
    gids.partition,
    gids.exchanger,
    gid_to_part,
    gids.ghost)
end

function PartitionedArrays.touched_hids(
  a::PRange,
  gids::AbstractPData{<:AbstractVector{<:Integer}})

  if a.gid_to_part !== nothing
    add_gids!(a,gids)
  end
  map_parts(touched_hids,a.partition,gids)
end

function remove_ghost_cells(model::DiscreteModel,gids::IndexSet)
  DiscreteModelPortion(model,gids.oid_to_lid)
end

function _exchange!(
  i_to_a::AbstractPData{<:Vector{<:Vector}},
  exchanger::Exchanger)

  n = length(get_part(i_to_a))
  for i in 1:n
    a = map_parts(i_to_a) do i_to_a
      i_to_a[i]
    end
    exchange!(a,exchanger)
  end
end

function _cut(bgmodel::DistributedDiscreteModel,args...)
  _cut(LevelSetCutter(),bgmodel,args...)
end

function _cut(cutter::Cutter,bgmodel::DistributedDiscreteModel,args...)
  gids = get_cell_gids(bgmodel)
  cuts = map_parts(local_views(bgmodel),local_views(gids)) do bgmodel,gids
    ownmodel = remove_ghost_cells(bgmodel,gids)
    cutgeo = cut(cutter,ownmodel,args...)
    change_bgmodel(cutgeo,bgmodel,gids.oid_to_lid)
  end
  ls_to_bgcell_to_inoutcut = map_parts(c->c.ls_to_bgcell_to_inoutcut,cuts)
  _exchange!(ls_to_bgcell_to_inoutcut,gids.exchanger)
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

function _aggregate(
  stragegy,
  cutgeo::DistributedEmbeddedDiscretization,
  in_or_out::Integer=IN)

  geo = get_part(map_parts(c->c.geo,local_views(cutgeo)))
  _aggregate(stragegy,cutgeo,geo,in_or_out)
end

function _aggregate(
  stragegy,
  cutgeo::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out::Integer=IN)

  bgf_to_ioc = compute_bgfacet_to_inoutcut(bgmodel,geo)
  _aggregate(stragegy,cutgeo,geo,in_or_out,bgf_to_ioc)
end

function _aggregate(
  stragegy,
  cutgeo::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out::Integer,
  bgf_to_ioc::AbstractPData{<:AbstractVector})

  map_parts(local_views(cutgeo),bgf_to_ioc) do cutgeo,bgf_to_ioc
    aggregate(stragegy,cutgeo,geo,in_or_out,bgf_to_ioc)
  end
  # TODO: parallel aggregation here
end

function change_bgmodel(
  cut::EmbeddedDiscretization,
  newmodel::DiscreteModel,
  cell_to_newcell)

  ls_to_bgc_to_ioc = map(cut.ls_to_bgcell_to_inoutcut) do bgc_to_ioc
    new_bgc_to_ioc = Vector{Int8}(undef,num_cells(newmodel))
    new_bgc_to_ioc[cell_to_newcell] = bgc_to_ioc
    new_bgc_to_ioc
  end
  subcells = change_bgmodel(cut.subcells,cell_to_newcell)
  subfacets = change_bgmodel(cut.subfacets,cell_to_newcell)
  EmbeddedDiscretization(
    newmodel,
    ls_to_bgc_to_ioc,
    subcells,
    cut.ls_to_subcell_to_inout,
    subfacets,
    cut.ls_to_subfacet_to_inout,
    cut.oid_to_ls,
    cut.geo)
end

function change_bgmodel(cells::SubCellData,cell_to_newcell)
  cell_to_bgcell = lazy_map(Reindex(cell_to_newcell),cells.cell_to_bgcell)
  SubCellData(
    cells.cell_to_points,
    collect(Int32,cell_to_bgcell),
    cells.point_to_coords,
    cells.point_to_rcoords)
end

function change_bgmodel(facets::SubFacetData,cell_to_newcell)
  facet_to_bgcell = lazy_map(Reindex(cell_to_newcell),facets.facet_to_bgcell)
  SubFacetData(
    facets.facet_to_points,
    facets.facet_to_normal,
    collect(Int32,facet_to_bgcell),
    facets.point_to_coords,
    facets.point_to_rcoords)
end

function compute_bgfacet_to_inoutcut(
  bgmodel::DistributedDiscreteModel,
  bgf_to_ioc::AbstractPData{<:AbstractVector})

  D = num_dims(get_part(local_views(bgmodel)))
  gids = get_cell_gids(bgmodel)
  bgf_to_ioc = map_parts(
    local_views(bgmodel),
    local_views(gids),
    bgf_to_ioc) do bgmodel,gids,bgf_to_ioc

    ownmodel = remove_ghost_cells(bgmodel,gids)
    f_to_pf = Gridap.Geometry.get_face_to_parent_face(ownmodel,D-1)
    _bgf_to_ioc = Vector{eltype(bgf_to_ioc)}(undef,num_faces(bgmodel,D-1))
    _bgf_to_ioc[f_to_pf] .= bgf_to_ioc
    _bgf_to_ioc
  end
  facet_gids = get_face_gids(bgmodel,D-1)
  exchange!(bgf_to_ioc,facet_gids.exchanger)
  bgf_to_ioc
end


function compute_bgfacet_to_inoutcut(
  cutter::Cutter,
  bgmodel::DistributedDiscreteModel,
  geo)

  gids = get_cell_gids(bgmodel)
  bgf_to_ioc = map_parts(local_views(bgmodel),local_views(gids)) do model,gids
    ownmodel = remove_ghost_cells(model,gids)
    compute_bgfacet_to_inoutcut(cutter,ownmodel,geo)
  end
  compute_bgfacet_to_inoutcut(bgmodel,bgf_to_ioc)
end


function compute_bgfacet_to_inoutcut(bgmodel::DistributedDiscreteModel,args...)
  cutter = LevelSetCutter()
  compute_bgfacet_to_inoutcut(cutter,bgmodel,args...)
end


# Driver

partition = (2,2)

parts = get_part_ids(SequentialBackend(),partition)

u(x) = x[1] - x[2]
f(x) = -Δ(u)(x)
ud(x) = u(x)

L = 1
p0 = Point(0.0,0.0)
pmin = p0-L/2
pmax = p0+L/2


R = 0.35
geo = disk(R,x0=p0)

R = 0.15
d = L/4
geo1 = disk(R,x0=p0-d)
geo2 = disk(R,x0=p0+d)
#geo = !union(geo1,geo2)

n = 8
mesh_partition = (n,n)
bgmodel = CartesianDiscreteModel(parts,pmin,pmax,mesh_partition)
bgtrian = Triangulation(bgmodel)
writevtk(bgtrian,"bgtrian")

dp = pmax - pmin
h = dp[1]/n

cutgeo = cut(bgmodel,geo)

strategy = AggregateAllCutCells()
strategy = AggregateCutCellsByThreshold(0.5)
aggregates = aggregate(strategy,cutgeo)

Ω_bg = Triangulation(bgmodel)
Ω_act = Triangulation(cutgeo,ACTIVE)
Ω = Triangulation(cutgeo,PHYSICAL)
Γ = EmbeddedBoundary(cutgeo)


writevtk(Ω_bg,"trian")
writevtk(Ω_act,"trian_act")
writevtk(Ω,"trian_O")
writevtk(Γ,"trian_G")

Γ = remove_ghost_cells(Γ)
Ω = remove_ghost_cells(Ω)

n_Γ = get_normal_vector(Γ)

order = 1
degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)

reffe = ReferenceFE(lagrangian,Float64,order)

Vstd = FESpace(Ω_act,reffe)

V = AgFEMSpace(bgmodel,Vstd,aggregates)
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
    "gid"=>_gids])#,
#  cellfields=["uh"=>uh])

writevtk(Ω,"trian_O",cellfields=["uh"=>uh])
writevtk(Γ,"trian_G")
@test el2/ul2 < 1.e-8
@test eh1/uh1 < 1.e-7


##

partition = (2,2)

parts = get_part_ids(SequentialBackend(),partition)

L = 1
p0 = Point(0.0,0.0)
pmin = p0-L/2
pmax = p0+L/2


R = 0.35
geo = disk(R,x0=p0)

n = 8
mesh_partition = (n,n)
bgmodel = CartesianDiscreteModel(parts,pmin,pmax,mesh_partition)

cutgeo = _cut(bgmodel,geo)

bgf_to_ioc = compute_bgfacet_to_inoutcut(bgmodel,geo)

# BUG: aggregates do not have information of ghost cut cells
# TODO: exchange measures and aggregates
_aggregate(strategy,cutgeo)


# TODO: parallel aggregation (parallel aggfem article)
# 0. exchange measures and aggregates through gid
# 1. root in ghost layer
# 2. root in neighbors
# 3. root in neighbors of neighbors

end # module
