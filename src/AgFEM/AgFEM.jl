module AgFEM

using LightGraphs
using LinearAlgebra

using Gridap
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces

using GridapEmbedded.CSG
using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: compute_bgcell_to_inoutcut

export aggregate
export color_aggregates
export AggregateAllCutCells
export AgFEMSpace

include("CellAggregation.jl")

include("AgFEMSpaces.jl")

end # module
