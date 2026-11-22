module LevelSetCutters

using AbstractTrees
import GridapEmbedded.CSG
using GridapEmbedded.CSG: Node, Leaf, replace_data
import GridapEmbedded.CSG: get_tree
import GridapEmbedded.CSG: similar_geometry
import GridapEmbedded.CSG: compatible_geometries

using GridapEmbedded.Interfaces
import GridapEmbedded.Interfaces: cut
import GridapEmbedded.Interfaces: cut_facets
import GridapEmbedded.Interfaces: compute_bgcell_to_inoutcut
import GridapEmbedded.Interfaces: compute_bgfacet_to_inoutcut
using GridapEmbedded.Interfaces: Simplex
using GridapEmbedded.Interfaces: UNDEF
using GridapEmbedded.Interfaces: merge_sub_face_data
using GridapEmbedded.Interfaces: compute_inoutcut
using GridapEmbedded.Interfaces: SubCellTriangulation, SubFacetTriangulation
using GridapEmbedded.Interfaces: SubFaceTriangulation, SubFacetBoundaryTriangulation

using LinearAlgebra
using MiniQhull
using StaticArrays
using FillArrays: Fill

using Gridap.TensorValues
using Gridap.ReferenceFEs
import Gridap.ReferenceFEs: num_vertices, num_facets, num_nodes, num_faces
import Gridap.ReferenceFEs: get_vertex_coordinates, get_face_vertices
using Gridap.Arrays
using Gridap.Arrays: testitem, return_type
using Gridap.Fields
using Gridap.Helpers
using Gridap.Geometry
using Gridap.Geometry: push_normal
using Gridap.CellData
using Gridap.Polynomials
using Gridap.Visualization
using Gridap.FESpaces
using Gridap.MultiField

using GridapDistributed: DistributedMultiFieldFESpace

export LevelSetCutter
export RecursivePlaneCutter
export AnalyticalGeometry
export doughnut
export popcorn
export tube
export olympic_rings
export sphere
export disk
export cylinder
export plane
export square
export quadrilateral
export cube
export discretize
export DiscreteGeometry

include("AnalyticalGeometries.jl")

include("DiscreteGeometries.jl")

include("LookupTables.jl")

include("CutTriangulations.jl")

include("DifferentiableTriangulations.jl")

# The recursive cell cutter and the cutter built on it
include("CellCutters.jl")

include("CutFaces.jl")

include("RecursivePlaneCutter.jl")

"""
    struct LevelSetCutter <: Cutter end

Cutter for `DiscreteGeometry` and `AnalyticalGeometry`.

## Usage

    cut(background::Union{DiscreteModel,Triangulation},geom::AnalyticalGeometry)
    cut(background::Union{DiscreteModel,Triangulation},geom::DiscreteGeometry)
    cut_facets(background::Union{DiscreteModel,Triangulation},geom::AnalyticalGeometry)
    cut_facets(background::Union{DiscreteModel,Triangulation},geom::DiscreteGeometry)

"""
struct LevelSetCutter <: Cutter end

# Default cutter dispatch
for f in (:cut,:cut_facets,:compute_bgcell_to_inoutcut,:compute_bgfacet_to_inoutcut)
  @eval $f(background::Union{DiscreteModel,Triangulation},geom::CSG.Geometry) =
    $f(LevelSetCutter(),background,geom)
end

# Unified DiscreteModel methods for Cutter subtypes
function cut(cutter::Cutter,background::DiscreteModel,geom::CSG.Geometry)
  data = _cut_grid(cutter,get_grid(background),geom,Val(true))
  EmbeddedDiscretization(background, data..., geom)
end

function cut_facets(cutter::Cutter,background::DiscreteModel,geom::CSG.Geometry)
  D = num_cell_dims(background)
  data = _cut_grid(cutter,Grid(ReferenceFE{D-1},background),geom,Val(false))
  EmbeddedFacetDiscretization(background, data..., geom)
end

function compute_bgcell_to_inoutcut(cutter::Cutter,model::DiscreteModel,geom::CSG.Geometry)
  _compute_bgcell_to_inoutcut(cutter,get_grid(model),geom)
end

function compute_bgfacet_to_inoutcut(cutter::Cutter,model::DiscreteModel,geom::CSG.Geometry)
  D = num_cell_dims(model)
  _compute_bgcell_to_inoutcut(cutter,Grid(ReferenceFE{D-1},model),geom)
end

# Unified subdomain Triangulation methods for Cutter subtypes
function cut(cutter::Cutter,trian::Triangulation,geom::CSG.Geometry)
  bgmodel, tcell_to_bgcell = _subdomain_cell_glue(trian)
  _is_full_model(tcell_to_bgcell,bgmodel) && return cut(cutter,bgmodel,geom)
  ls_to_ioc, subcells, ls_to_cio, subfacets, ls_to_fio, oid_to_ls =
    _cut_grid(cutter,_get_subdomain_grid(trian),geom,Val(true))
  n_bgcells = num_cells(bgmodel)
  ls_to_bgioc = map(ioc -> _extend_inoutcut(ioc,tcell_to_bgcell,n_bgcells),ls_to_ioc)
  EmbeddedDiscretization(
    bgmodel,ls_to_bgioc,
    _remap_bgcell(subcells,tcell_to_bgcell),ls_to_cio,
    _remap_bgcell(subfacets,tcell_to_bgcell),ls_to_fio,
    oid_to_ls,geom)
end

function cut_facets(cutter::Cutter,trian::Triangulation,geom::CSG.Geometry)
  bgmodel, tcell_to_bgcell = _subdomain_cell_glue(trian)
  _is_full_model(tcell_to_bgcell,bgmodel) && return cut_facets(cutter,bgmodel,geom)
  tfacet_grid, tfacet_to_bgfacet, n_bgfacets = _subdomain_facet_grid(trian,bgmodel,tcell_to_bgcell)
  ls_to_ioc, subfacets, ls_to_fio, oid_to_ls = _cut_grid(cutter,tfacet_grid,geom,Val(false))
  ls_to_bgioc = map(ioc -> _extend_inoutcut(ioc,tfacet_to_bgfacet,n_bgfacets),ls_to_ioc)
  EmbeddedFacetDiscretization(
    bgmodel,ls_to_bgioc,
    _remap_bgcell(subfacets,tfacet_to_bgfacet),ls_to_fio,
    oid_to_ls,geom)
end

function compute_bgcell_to_inoutcut(cutter::Cutter,trian::Triangulation,geom::CSG.Geometry)
  bgmodel, tcell_to_bgcell = _subdomain_cell_glue(trian)
  _is_full_model(tcell_to_bgcell,bgmodel) && return compute_bgcell_to_inoutcut(cutter,bgmodel,geom)
  tcell_to_ioc = _compute_bgcell_to_inoutcut(cutter,_get_subdomain_grid(trian),geom)
  _extend_inoutcut(tcell_to_ioc,tcell_to_bgcell,num_cells(bgmodel))
end

function compute_bgfacet_to_inoutcut(cutter::Cutter,trian::Triangulation,geom::CSG.Geometry)
  bgmodel, tcell_to_bgcell = _subdomain_cell_glue(trian)
  _is_full_model(tcell_to_bgcell,bgmodel) && return compute_bgfacet_to_inoutcut(cutter,bgmodel,geom)
  tfacet_grid, tfacet_to_bgfacet, n_bgfacets = _subdomain_facet_grid(trian,bgmodel,tcell_to_bgcell)
  tfacet_to_ioc = _compute_bgcell_to_inoutcut(cutter,tfacet_grid,geom)
  _extend_inoutcut(tfacet_to_ioc,tfacet_to_bgfacet,n_bgfacets)
end

function _cut_domain(cut::EmbeddedDiscretization)
  if !isempty(cut.ls_to_bgcell_to_inoutcut)
    mask = map(!=(UNDEF),first(cut.ls_to_bgcell_to_inoutcut))
    !all(mask) && return Triangulation(cut.bgmodel,mask)
  end
  cut.bgmodel
end

cut_facets(cut::EmbeddedDiscretization,geo::CSG.Geometry) = cut_facets(_cut_domain(cut),geo)
cut_facets(cutter::Cutter,cut::EmbeddedDiscretization,geo::CSG.Geometry) = cut_facets(cutter,_cut_domain(cut),geo)
cut_facets(cut::EmbeddedFacetDiscretization,args...) = cut

# LevelSetCutter grid-level hooks
function _cut_grid(::LevelSetCutter,grid::Grid,geom,::Val{true})
  subcells0, ls_to_point_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls = initial_sub_triangulation(grid,geom)
  subcells, ls_to_cell_to_inout, subfacets, ls_to_facet_to_inout =
    cut_sub_triangulation_with_boundary_several_levelsets(subcells0,ls_to_point_to_value)
  ls_to_bgcell_to_inoutcut, subcells, ls_to_cell_to_inout, subfacets, ls_to_facet_to_inout, oid_to_ls
end

function _cut_grid(::LevelSetCutter,grid::Grid,geom,::Val{false})
  subcells0, ls_to_point_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls = initial_sub_triangulation(grid,geom)
  subcells, ls_to_cell_to_inout = cut_sub_triangulation_several_levelsets(subcells0,ls_to_point_to_value)
  ls_to_bgcell_to_inoutcut, subcells, ls_to_cell_to_inout, oid_to_ls
end

_cut_grid(cutter::RecursivePlaneCutter,grid::Grid,geom,wb::Val) = _cut(cutter,grid,geom,wb)

# Shared subdomain helpers
function _subdomain_cell_glue(trian::Triangulation)
  bgmodel = get_background_model(trian)
  D = num_cell_dims(bgmodel)
  @check num_cell_dims(trian) == D """
  Subdomain cutting requires a Triangulation of the same cell dimension ($D) as its background DiscreteModel.
  """
  glue = get_glue(trian,Val(D))
  @check isa(glue,FaceToFaceGlue) && !isnothing(glue.mface_to_tface) """
  Subdomain cutting requires an injective body-fitted Triangulation of the background DiscreteModel.
  """
  return bgmodel, glue.tface_to_mface
end

_is_full_model(tcell_to_bgcell,bgmodel) =
  isa(tcell_to_bgcell,IdentityVector) && length(tcell_to_bgcell) == num_cells(bgmodel)

function _extend_inoutcut(local_ioc::AbstractVector{Int8},local_to_global,n_global::Integer)
  global_ioc = fill(Int8(UNDEF),n_global)
  global_ioc[local_to_global] .= local_ioc
  return global_ioc
end

function _remap_bgcell(st::SubFaceData{Df},local_to_global) where Df
  cell_to_bgcell = collect(Int32,lazy_map(Reindex(local_to_global),st.cell_to_bgcell))
  SubFaceData{Df}(st.cell_to_points,cell_to_bgcell,st.point_to_coords,st.point_to_rcoords,st.cell_to_normal)
end

function _subdomain_facet_grid(trian::Triangulation,bgmodel::DiscreteModel,tcell_to_bgcell)
  D = num_cell_dims(bgmodel)
  cell_to_facets = Table(get_faces(get_grid_topology(bgmodel),D,D-1))
  n_bgfacets = num_facets(bgmodel)
  bgfacet_in_trian = fill(false,n_bgfacets)
  for bgcell in tcell_to_bgcell
    for p in cell_to_facets.ptrs[bgcell]:(cell_to_facets.ptrs[bgcell+1]-1)
      bgfacet_in_trian[cell_to_facets.data[p]] = true
    end
  end
  tfacet_to_bgfacet = findall(bgfacet_in_trian)
  tfacet_grid = GridPortion(Grid(ReferenceFE{D-1},bgmodel),tfacet_to_bgfacet)
  cgrid = _get_subdomain_grid(trian)
  if isa(cgrid,Geometry.GridPortion)
    @check tfacet_grid.node_to_parent_node == cgrid.node_to_parent_node
  end
  return tfacet_grid, tfacet_to_bgfacet, n_bgfacets
end

function _compute_bgcell_to_inoutcut(::LevelSetCutter,grid::Grid,geom)
  ls_to_cell_to_inoutcut, tree, oid_to_ls = _compute_ls_to_bgcell_to_inoutcut(grid,geom)

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    ls = oid_to_ls[oid]
    cell_to_inoutcut = ls_to_cell_to_inoutcut[ls]
    cell_to_inoutcut, name, meta
  end

  newtree = replace_data(identity,conversion,tree)
  compute_inoutcut(newtree)
end

function _compute_ls_to_bgcell_to_inoutcut(grid,geom::AnalyticalGeometry)
  _compute_ls_to_bgcell_to_inoutcut(grid,discretize(geom,grid))
end

function _compute_ls_to_bgcell_to_inoutcut(grid,geom::DiscreteGeometry)
  tree = get_tree(geom)
  ls_to_point_to_value, oid_to_ls = _find_unique_leaves(tree)
  ls_to_point_to_value = _adapt_point_values_to_grid(ls_to_point_to_value,grid)
  cell_to_points = get_cell_node_ids(UnstructuredGrid(grid))
  ls_to_cell_to_inoutcut = [
    _compute_in_out_or_cut(cell_to_points,point_to_value)
    for point_to_value in ls_to_point_to_value]
  ls_to_cell_to_inoutcut, tree, oid_to_ls
end

end #module
