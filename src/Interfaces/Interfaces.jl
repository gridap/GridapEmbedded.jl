module Interfaces

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Visualization

import Gridap.Geometry: UnstructuredGrid
import Gridap.Visualization: writevtk

export IN
export OUT
export INTERFACE
export CUT
export EmbeddedDiscretization
export SubTriangulation
export FacetSubTriangulation
export Cutter
export cut
export split_in_out

include("EmbeddedDiscretizations.jl")

include("Cutters.jl")

end # module
