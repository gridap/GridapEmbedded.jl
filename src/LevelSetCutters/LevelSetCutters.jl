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
using GridapEmbedded.Interfaces: merge_sub_face_data
using GridapEmbedded.Interfaces: compute_inoutcut

using LinearAlgebra
using MiniQhull

using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Fields
using Gridap.Helpers
using Gridap.Geometry
using Gridap.CellData
using Gridap.Polynomials
using Gridap.Visualization

export LevelSetCutter
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

struct LevelSetCutter <: Cutter end

function cut(cutter::LevelSetCutter,background::DiscreteModel,geom)
  data = _cut_ls(background,geom)
  EmbeddedDiscretization(background, data..., geom)
end

function cut(background::DiscreteModel,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  cut(cutter,background,geom)
end

function cut(background::DiscreteModel,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  cut(cutter,background,geom)
end

function _cut_ls(model::DiscreteModel,geom)
  grid = get_grid(model)
  _cut_ls(grid,geom)
end

function _cut_ls(grid::Grid,geom)
  out = initial_sub_triangulation(grid,geom)
  subcells0, ls_to_point_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls = out
  out2 = cut_sub_triangulation_with_boundary_several_levelsets(subcells0,ls_to_point_to_value)
  subcells, ls_to_cell_to_inout, subfacets, ls_to_facet_to_inout  = out2
  ls_to_bgcell_to_inoutcut, subcells, ls_to_cell_to_inout, subfacets, ls_to_facet_to_inout, oid_to_ls
end

function cut_facets(cutter::LevelSetCutter,background::DiscreteModel,geom)
  data = _cut_ls_facets(background,geom)
  EmbeddedFacetDiscretization(background, data..., geom)
end

function cut_facets(background::DiscreteModel,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  cut_facets(cutter,background,geom)
end

function cut_facets(background::DiscreteModel,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  cut_facets(cutter,background,geom)
end

function _cut_ls_facets(model::DiscreteModel,geom)
  D = num_cell_dims(model)
  grid = Grid(ReferenceFE{D-1},model)
  _cut_ls_facets(grid,geom)
end

function _cut_ls_facets(grid::Grid,geom)
  out = initial_sub_triangulation(grid,geom)
  subcells0, ls_to_point_to_value, ls_to_bgcell_to_inoutcut, oid_to_ls = out
  out2 = cut_sub_triangulation_several_levelsets(subcells0,ls_to_point_to_value)
  subcells, ls_to_cell_to_inout  = out2
  ls_to_bgcell_to_inoutcut, subcells, ls_to_cell_to_inout, oid_to_ls
end

function compute_bgcell_to_inoutcut(model::DiscreteModel,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  compute_bgcell_to_inoutcut(cutter,model,geom)
end

function compute_bgcell_to_inoutcut(model::DiscreteModel,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  compute_bgcell_to_inoutcut(cutter,model,geom)
end

function compute_bgcell_to_inoutcut(::LevelSetCutter,model::DiscreteModel,geom)
  grid = get_grid(model)
  _compute_bgcell_to_inoutcut(grid,geom)
end

function compute_bgfacet_to_inoutcut(model::DiscreteModel,geom::AnalyticalGeometry)
  cutter = LevelSetCutter()
  compute_bgfacet_to_inoutcut(cutter,model,geom)
end

function compute_bgfacet_to_inoutcut(model::DiscreteModel,geom::DiscreteGeometry)
  cutter = LevelSetCutter()
  compute_bgfacet_to_inoutcut(cutter,model,geom)
end

function compute_bgfacet_to_inoutcut(::LevelSetCutter,model::DiscreteModel,geom)
  D = num_cell_dims(model)
  grid = Grid(ReferenceFE{D-1},model)
  _compute_bgcell_to_inoutcut(grid,geom)
end

function _compute_bgcell_to_inoutcut(grid,geom)
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
  ugrid = UnstructuredGrid(grid)
  tree = get_tree(geom)
  ls_to_point_to_value, oid_to_ls = _find_unique_leaves(tree)
  p = _check_and_get_polytope(ugrid)
  table = LookupTable(p)
  cell_to_points = get_cell_node_ids(ugrid)
  ls_to_cell_to_inoutcut = [
    _compute_in_out_or_cut(table,cell_to_points,point_to_value)
    for point_to_value in ls_to_point_to_value]
  ls_to_cell_to_inoutcut, tree, oid_to_ls
end

end #module
