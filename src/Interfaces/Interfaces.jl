module Interfaces

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.Visualization

import Gridap.Geometry: UnstructuredGrid
import Gridap.Geometry: BoundaryTriangulation
import Gridap.Geometry: SkeletonTriangulation
import Gridap.Geometry: Triangulation
import Gridap.Geometry: DiscreteModel
import Gridap.Geometry: get_node_coordinates
import Gridap.Geometry: get_cell_node_ids
import Gridap.Geometry: get_cell_coordinates
import Gridap.Geometry: get_reffes
import Gridap.Geometry: get_cell_type
import Gridap.Geometry: get_cell_to_bgcell
import Gridap.Geometry: TriangulationStyle
import Gridap.Geometry: get_background_triangulation
import Gridap.Geometry: get_cell_ref_map
import Gridap.Geometry: get_facet_normal
import Gridap.Geometry: compress_contributions
import Gridap.Geometry: compress_ids

using GridapEmbedded.CSG

export IN
export OUT
export INTERFACE
export CUT
export CUTIN
export CUTOUT
export EmbeddedDiscretization
export EmbeddedFacetDiscretization
export SubCellData
export SubFacetData
export Cutter
export cut
export cut_facets
export compute_bgcell_to_inoutcut
export compute_bgfacet_to_inoutcut
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

include("SubCellTriangulations.jl")

include("SubFacetTriangulations.jl")

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
