module GridapEmbedded

using MiniQhull
import MiniQhull: delaunay
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Fields
using Gridap.Helpers
using Gridap.Geometry
using Gridap.Integration

export LookupTable
export num_cases
export compute_case

include("LookupTables.jl")

end # module
