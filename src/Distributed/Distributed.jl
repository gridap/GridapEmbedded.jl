module Distributed

using Gridap
using GridapDistributed
using PartitionedArrays
using FillArrays

using Gridap.Arrays
using Gridap.CellData
using Gridap.Geometry
using Gridap.Helpers
using Gridap.ReferenceFEs
using Gridap.FESpaces

using PartitionedArrays: VectorFromDict

using GridapEmbedded.CSG
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: Cutter
using GridapEmbedded.Interfaces: ActiveInOrOut
using GridapEmbedded.Interfaces: SubFacetTriangulation
using GridapEmbedded.Interfaces: SubCellData
using GridapEmbedded.Interfaces: SubFacetData
using GridapEmbedded.Interfaces: AbstractEmbeddedDiscretization
using GridapEmbedded.AgFEM: _touch_aggregated_cells!
using GridapEmbedded.AgFEM: AggregateCutCellsByThreshold
using GridapEmbedded.MomentFittedQuadratures: MomentFitted
using Gridap.Geometry: AppendedTriangulation, TriangulationView
using Gridap.Geometry: get_face_to_parent_face
using Gridap.Arrays: find_inverse_index_map, testitem, return_type
using Gridap.FESpaces: FESpaceWithLinearConstraints
using Gridap.FESpaces: _dof_to_DOF, _DOF_to_dof

using GridapDistributed: DistributedDiscreteModel, DistributedTriangulation, DistributedMeasure
using GridapDistributed: DistributedFESpace, DistributedSingleFieldFESpace
using GridapDistributed: NoGhost, WithGhost, filter_cells_when_needed, add_ghost_cells
using GridapDistributed: generate_gids, generate_cell_gids
using GridapDistributed: _find_vector_type

import GridapEmbedded.AgFEM: aggregate
import GridapEmbedded.AgFEM: AgFEMSpace
import GridapEmbedded.Interfaces: cut
import GridapEmbedded.Interfaces: cut_facets
import GridapEmbedded.Interfaces: EmbeddedBoundary
import GridapEmbedded.Interfaces: compute_bgfacet_to_inoutcut
import GridapEmbedded.Interfaces: compute_bgcell_to_inoutcut
import GridapEmbedded.Interfaces: GhostSkeleton
import GridapEmbedded.CSG: get_geometry
import GridapEmbedded.LevelSetCutters: discretize, DiscreteGeometry
import Gridap.Geometry: Triangulation
import Gridap.Geometry: SkeletonTriangulation
import Gridap.Geometry: BoundaryTriangulation
import Gridap.Geometry: get_background_model
import GridapDistributed: local_views
import GridapDistributed: remove_ghost_cells

include("DistributedDiscretizations.jl")

include("DistributedDiscreteGeometries.jl")

include("DistributedSubFacetTriangulations.jl")

include("DistributedAgFEM.jl")

include("DistributedQuadratures.jl")

end # module
