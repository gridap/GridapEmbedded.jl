module AgFEM

using Graphs
using LinearAlgebra

using Gridap
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Geometry: get_cell_to_parent_cell
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.FESpaces

using GridapEmbedded.CSG
using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: compute_subcell_to_inout

export aggregate
export color_aggregates
export aggregate_narrow_band
export AggregateAllCutCells
export compute_cell_bboxes
export compute_cell_to_dface_bboxes
export AgFEMSpace

include("CellAggregation.jl")

include("AggregateBoundingBoxes.jl")

include("AgFEMSpaces.jl")

end # module
