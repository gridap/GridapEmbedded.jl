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
export aggregatespace
export color_aggregates
export color_aggregates_space
export AggregateAllCutCells
export AggregateSpaceCutCells
export AgFEMSpace

include("CellAggregation.jl")

include("AgFEMSpaces.jl")

include("CellAggregation_space.jl")
end # module
