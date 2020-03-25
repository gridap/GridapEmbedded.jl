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

include("LookupTables.jl")

include("Geometries.jl")

include("SubTriangulations.jl")

struct LevelSetCutter <: Cutter end

function cut(cutter::LevelSetCutter,background,geom)
  data = _cut_ls(background,geom)
  EmbeddedDiscretization(background, data...)
end

function _cut_ls(model::DiscreteModel,geom)
  grid = get_grid(model)
  _cut_ls(grid,geom)
end

function _cut_ls(grid::Grid,geom)
  subtrian, subgeom, bgcell_to_inoutcut = initial_sub_triangulation(grid,geom)
  st, ls_to_fst = cut_sub_triangulation(subtrian,subgeom)
  ls_to_name = [ "tag_$i" for i in 1:length(ls_to_fst)]
  bgcell_to_inoutcut, st, ls_to_fst, ls_to_name
end

end #module


