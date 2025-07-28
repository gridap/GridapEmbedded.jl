module Interfaces

using FillArrays

using Gridap
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.Visualization

import GridapEmbedded.CSG: get_geometry
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
import Gridap.Geometry: get_background_model
import Gridap.Geometry: get_active_model
import Gridap.Geometry: compute_active_model
import Gridap.Geometry: get_glue
import Gridap.Geometry: get_grid
import Gridap.Geometry: FaceToFaceGlue
import Gridap.Geometry: get_facet_normal
import Gridap.Geometry: move_contributions
import Gridap.Geometry: is_change_possible
using Gridap.Geometry: GenericTriangulation
using Gridap.Geometry: CompositeTriangulation
using Gridap.Geometry: TriangulationView
using Gridap.Geometry: restrict

import Gridap.CellData: get_normal_vector
import Gridap.CellData: get_tangent_vector
import Gridap.Geometry: get_facet_normal

using GridapEmbedded.CSG

export IN
export OUT
export INTERFACE
export CUT
export CUT_IN
export CUT_OUT
export ACTIVE
export ACTIVE_IN
export ACTIVE_OUT
export PHYSICAL
export PHYSICAL_IN
export PHYSICAL_OUT
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
const CUT_IN = CutInOrOut(IN)
const CUT_OUT = CutInOrOut(OUT)
const PHYSICAL_IN = (CUT_IN,IN)
const PHYSICAL_OUT = (CUT_OUT,OUT)
const PHYSICAL = PHYSICAL_IN

struct ActiveInOrOut
  in_or_out::Int
end
const ACTIVE_IN = ActiveInOrOut(IN)
const ACTIVE_OUT = ActiveInOrOut(OUT)
const ACTIVE = ACTIVE_IN

include("SubCellTriangulations.jl")

include("SubFacetTriangulations.jl")

include("EmbeddedDiscretizations.jl")

include("EmbeddedFacetDiscretizations.jl")

include("CutFaceBoundaryTriangulations.jl")

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
