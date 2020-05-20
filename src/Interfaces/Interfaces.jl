module Interfaces

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Visualization

import Gridap.Geometry: UnstructuredGrid
import Gridap.Geometry: BoundaryTriangulation
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

using GridapEmbedded.CSG

export IN
export OUT
export INTERFACE
export CUT
export CUTIN
export CUTOUT
export EmbeddedDiscretization
export EmbeddedFacetDiscretization
export SubTriangulation
export FacetSubTriangulation
export Cutter
export cut
export cut_facets
export split_in_out
export EmbeddedBoundary
export GhostSkeleton

const IN = -1
const OUT = 1
const INTERFACE = 0
const CUT = 0

struct CutInOrOut
  in_or_out::Int
end
const CUTIN = CutInOrOut(IN)
const CUTOUT = CutInOrOut(OUT)

include("SubTriangulations.jl")

include("FacetSubTriangulations.jl")

include("EmbeddedDiscretizations.jl")

include("EmbeddedFacetDiscretizations.jl")

include("Cutters.jl")

function Simplex(p::Polytope)
  D = num_cell_dims(p)
  Simplex(Val{D}())
end

function Simplex(::Val{D}) where D
  extrusion = tfill(TET_AXIS,Val{D}())
  ExtrusionPolytope(extrusion)
end

function Simplex(::Val{2})
  TRI
end

function Simplex(::Val{3})
  TET
end

end # module
