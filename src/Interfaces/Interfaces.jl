module Interfaces

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Visualization

import Gridap.Geometry: UnstructuredGrid
import Gridap.Geometry: Triangulation
import Gridap.Visualization: writevtk
import Gridap.Visualization: DiscreteModel
import Gridap.Geometry: get_node_coordinates
import Gridap.Geometry: get_cell_nodes
import Gridap.Geometry: get_cell_coordinates
import Gridap.Geometry: get_reffes
import Gridap.Geometry: get_cell_type
import Gridap.Geometry: get_normal_vector
import Gridap.Geometry: get_face_to_cell
import Gridap.Geometry: get_face_to_cell_map
import Gridap.Geometry: restrict
import Gridap.Geometry: get_cell_id

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
export EmbeddedBoundary
export GhostSkeleton
export cell_measure

const IN = -1
const OUT = 1
const INTERFACE = 0
const CUT = 0

include("SubTriangulations.jl")

include("FacetSubTriangulations.jl")

include("EmbeddedDiscretizations.jl")

include("Cutters.jl")

end # module
