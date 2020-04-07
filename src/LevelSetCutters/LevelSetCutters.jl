module LevelSetCutters

using AbstractTrees
import GridapEmbedded.CSG
using GridapEmbedded.CSG: Node, Leaf, replace_data
import GridapEmbedded.CSG: get_tree
import GridapEmbedded.CSG: similar_geometry
import GridapEmbedded.CSG: compatible_geometries

using GridapEmbedded.Interfaces
import GridapEmbedded.Interfaces: cut
using GridapEmbedded.Interfaces: Simplex

using LinearAlgebra
using MiniQhull

using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Fields
using Gridap.Helpers
using Gridap.Geometry
using Gridap.Integration
using Gridap.Polynomials
using Gridap.Visualization

export LevelSetCutter
export doughnut
export tube
export olympic_rings
export sphere
export disk
export discretize

include("AnalyticalGeometries.jl")

include("DiscreteGeometries.jl")

include("LookupTables.jl")

include("SubTriangulations.jl")

struct LevelSetCutter <: Cutter end

function cut(cutter::LevelSetCutter,background::DiscreteModel,geom)
  data = _cut_ls(background,geom)
  EmbeddedDiscretization(background, data..., geom)
end

function cut(background::DiscreteModel,geom::AnalyticalGeometry)
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
  out2 = cut_sub_triangulation(subcells0,ls_to_point_to_value)
  subcells, ls_to_cell_to_inout, ls_to_subfacets, ls_to_ls_to_facet_to_inout  = out2
  ls_to_bgcell_to_inoutcut, subcells, ls_to_cell_to_inout, ls_to_subfacets, ls_to_ls_to_facet_to_inout, oid_to_ls
end

end #module



