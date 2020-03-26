module LevelSetCutters

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
export disc

include("LookupTables.jl")

include("Geometries.jl")

include("SubTriangulations.jl")

struct LevelSetCutter <: Cutter end

function cut(cutter::LevelSetCutter,background::DiscreteModel,geom)
  data = _cut_ls(background,geom)
  EmbeddedDiscretization(background, data...)
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
  subtrian, subgeom, bgcell_to_inoutcut = initial_sub_triangulation(grid,geom)
  st, ls_to_fst = cut_sub_triangulation(subtrian,subgeom)
  st_in, st_out = split_in_out(st)
  ls_to_name = reverse(geom.ls_to_name)
  bgcell_to_inoutcut, st_in, st_out, ls_to_fst, ls_to_name
end

end #module



