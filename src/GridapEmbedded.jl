module GridapEmbedded

using LinearAlgebra
using MiniQhull
import MiniQhull: delaunay
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

export LookupTable
export num_cases
export compute_case
export initial_sub_triangulation
export cut_sub_triangulation

include("LookupTables.jl")

include("SubTriangulations.jl")

end # module
