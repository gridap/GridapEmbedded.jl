module GridapEmbedded

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

import Gridap.Geometry: UnstructuredGrid
import Gridap.Visualization: writevtk

include("LookupTables.jl")

include("Geometries.jl")

include("SubTriangulations.jl")

end # module
