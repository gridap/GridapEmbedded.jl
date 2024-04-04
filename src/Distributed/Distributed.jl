module Distributed

using Gridap
using GridapDistributed
using PartitionedArrays

using Gridap.Arrays
using Gridap.Geometry
using Gridap.Helpers

using GridapEmbedded.CSG
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: Cutter
using GridapEmbedded.Interfaces: ActiveInOrOut
using GridapEmbedded.Interfaces: SubFacetTriangulation
using GridapEmbedded.Interfaces: SubCellData
using GridapEmbedded.Interfaces: SubFacetData
using Gridap.Geometry: AppendedTriangulation
using GridapDistributed: DistributedDiscreteModel
using GridapDistributed: DistributedTriangulation
using GridapDistributed: DistributedFESpace
using GridapDistributed: DistributedSingleFieldFESpace
using GridapDistributed: add_ghost_cells
using GridapDistributed: generate_gids
using GridapDistributed: generate_cell_gids
using GridapDistributed: _find_vector_type

import GridapEmbedded.AgFEM: aggregate
import GridapEmbedded.AgFEM: AgFEMSpace
import GridapEmbedded.Interfaces: cut
import GridapEmbedded.Interfaces: EmbeddedBoundary
import GridapEmbedded.Interfaces: compute_bgfacet_to_inoutcut
import GridapEmbedded.Interfaces: compute_bgcell_to_inoutcut
import Gridap.Geometry: Triangulation
import Gridap.Geometry: get_background_model
import GridapDistributed: local_views
import GridapDistributed: remove_ghost_cells


struct DistributedEmbeddedDiscretization{Dp,T,A,B} <: GridapType
  discretizations::A
  model::B
  function DistributedEmbeddedDiscretization(
    d::AbstractArray{<:EmbeddedDiscretization{Dp,T}},
    model::DistributedDiscreteModel) where {Dp,T}
    A = typeof(d)
    B = typeof(model)
    new{Dp,T,A,B}(d,model)
  end
end

local_views(a::DistributedEmbeddedDiscretization) = a.discretizations

get_background_model(a::DistributedEmbeddedDiscretization) = a.model

function cut(bgmodel::DistributedDiscreteModel,args...)
  cut(LevelSetCutter(),bgmodel,args...)
end

function cut(cutter::Cutter,bgmodel::DistributedDiscreteModel,args...)
  gids = get_cell_gids(bgmodel)
  cuts = map(local_views(bgmodel),local_views(gids)) do bgmodel,gids
    ownmodel = remove_ghost_cells(bgmodel,gids)
    cutgeo = cut(cutter,ownmodel,args...)
    change_bgmodel(cutgeo,bgmodel,own_to_local(gids))
  end
  ls_to_bgcell_to_inoutcut = map(c->c.ls_to_bgcell_to_inoutcut,cuts)
  _consistent!(ls_to_bgcell_to_inoutcut,gids)
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

function Triangulation(
  cutgeo::DistributedEmbeddedDiscretization,
  in_or_out::ActiveInOrOut,
  args...)

  distributed_embedded_triangulation(Triangulation,cutgeo,in_or_out,args...)
end

function Triangulation(cutgeo::DistributedEmbeddedDiscretization,args...)
  trian = distributed_embedded_triangulation(Triangulation,cutgeo,args...)
  remove_ghost_cells(trian)
end

function EmbeddedBoundary(cutgeo::DistributedEmbeddedDiscretization,args...)
  trian = distributed_embedded_triangulation(EmbeddedBoundary,cutgeo,args...)
  remove_ghost_cells(trian)
end

function distributed_embedded_triangulation(
  T,
  cutgeo::DistributedEmbeddedDiscretization,
  args...)

  trians = map(local_views(cutgeo)) do lcutgeo
    T(lcutgeo,args...)
  end
  bgmodel = get_background_model(cutgeo)
  DistributedTriangulation(trians,bgmodel)
end


function aggregate(stragegy,cutgeo::DistributedEmbeddedDiscretization)
  aggregates = map(local_views(cutgeo)) do lcutgeo
    aggregate(stragegy,lcutgeo)
  end
  bgmodel = get_background_model(cutgeo)
  consistent_aggregates(aggregates,bgmodel)
end

function consistent_aggregates(aggregates,bgmodel::DistributedDiscreteModel)
  consistent_aggregates(aggregates,get_cell_gids(bgmodel))
end

function consistent_aggregates(aggregates,gids::PRange)
  global_aggregates = map(aggregates,local_to_global(gids)) do agg,gid
    map(i-> iszero(i) ? i : gid[i],agg)
  end
  paggregates = PVector(global_aggregates,partition(gids))
  consistent!(paggregates) |> wait
  map(local_values(paggregates),global_to_local(gids)) do agg,lgid
    map(i-> iszero(i) ? i : lgid[i],agg)
  end
end

function AgFEMSpace(
  bgmodel::DistributedDiscreteModel,
  f::DistributedFESpace,
  bgcell_to_bgcellin::AbstractArray{<:AbstractVector},
  g::DistributedFESpace=f)

  bgmodel_gids = get_cell_gids(bgmodel)
  spaces = map(
    local_views(f),
    bgcell_to_bgcellin,
    local_views(g),
    local_views(bgmodel_gids)) do f,bgcell_to_bgcellin,g,gids
      AgFEMSpace(f,bgcell_to_bgcellin,g,local_to_global(gids))
  end
  trians = map(get_triangulation,local_views(f))
  trian = DistributedTriangulation(trians,bgmodel)
  trian = add_ghost_cells(trian)
  trian_gids = generate_cell_gids(trian)
  cell_to_ldofs = map(get_cell_dof_ids,spaces)
  cell_to_ldofs = map(i->map(sort,i),cell_to_ldofs)
  nldofs = map(num_free_dofs,spaces)
  gids = generate_gids(trian_gids,cell_to_ldofs,nldofs)
  vector_type = _find_vector_type(spaces,gids)
  DistributedSingleFieldFESpace(spaces,gids,vector_type)
end

function remove_ghost_cells(trian::DistributedTriangulation)
  model = get_background_model(trian)
  gids = get_cell_gids(model)
  trians = map(local_views(trian),local_views(gids)) do trian,gids
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

function remove_ghost_cells(model::DiscreteModel,gids::AbstractLocalIndices)
  DiscreteModelPortion(model,own_to_local(gids))
end

function _consistent!(
  p_to_i_to_a::AbstractArray{<:Vector{<:Vector}},
  prange::PRange)

  tables = map(Table,p_to_i_to_a)
  data = map(t->t.data,tables)
  pdata = PVector(data,partition(prange))
  consistent!(pdata) |> wait
  map(copyto!,p_to_i_to_a,tables)
end


function compute_bgfacet_to_inoutcut(
  bgmodel::DistributedDiscreteModel,
  bgf_to_ioc::AbstractArray{<:AbstractVector})

  D = num_dims(eltype(local_views(bgmodel)))
  gids = get_cell_gids(bgmodel)
  bgf_to_ioc = map(
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
  pbgf_to_ioc = PVector(bgf_to_ioc,partition(facet_gids))
  consistent!(pbgf_to_ioc) |> wait
  local_values(pbgf_to_ioc)
end


function compute_bgfacet_to_inoutcut(
  cutter::Cutter,
  bgmodel::DistributedDiscreteModel,
  geo)

  gids = get_cell_gids(bgmodel)
  bgf_to_ioc = map(local_views(bgmodel),local_views(gids)) do model,gids
    ownmodel = remove_ghost_cells(model,gids)
    compute_bgfacet_to_inoutcut(cutter,ownmodel,geo)
  end
  compute_bgfacet_to_inoutcut(bgmodel,bgf_to_ioc)
end


function compute_bgfacet_to_inoutcut(bgmodel::DistributedDiscreteModel,args...)
  cutter = LevelSetCutter()
  compute_bgfacet_to_inoutcut(cutter,bgmodel,args...)
end

function compute_bgcell_to_inoutcut(cutgeo::DistributedEmbeddedDiscretization,args...)
  map(local_views(cutgeo)) do cutgeo
    compute_bgcell_to_inoutcut(cutgeo,args...)
  end
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

end # module
