module Distributed

using Gridap
using GridapDistributed
using PartitionedArrays
using FillArrays

using Gridap.Arrays
using Gridap.Geometry
using Gridap.Helpers
using Gridap.ReferenceFEs

using GridapEmbedded.CSG
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: Cutter
using GridapEmbedded.Interfaces: ActiveInOrOut
using GridapEmbedded.Interfaces: SubFacetTriangulation
using GridapEmbedded.Interfaces: SubCellData
using GridapEmbedded.Interfaces: SubFacetData
using GridapEmbedded.AgFEM: _touch_aggregated_cells!
using GridapEmbedded.AgFEM: AggregateCutCellsByThreshold
using Gridap.Geometry: AppendedTriangulation
using GridapDistributed: DistributedDiscreteModel
using GridapDistributed: DistributedTriangulation
using GridapDistributed: DistributedFESpace
using GridapDistributed: DistributedSingleFieldFESpace
using GridapDistributed: add_ghost_cells
using GridapDistributed: generate_gids
using GridapDistributed: generate_cell_gids
using GridapDistributed: _find_vector_type

import GridapEmbedded.AgFEM: aggregate
import GridapEmbedded.AgFEM: AgFEMSpace
import GridapEmbedded.Interfaces: cut
import GridapEmbedded.Interfaces: EmbeddedBoundary
import GridapEmbedded.Interfaces: compute_bgfacet_to_inoutcut
import GridapEmbedded.Interfaces: compute_bgcell_to_inoutcut
import Gridap.Geometry: Triangulation
import Gridap.Geometry: get_background_model
import GridapDistributed: local_views
import GridapDistributed: remove_ghost_cells

include("DistributedDiscretizations.jl")

include("DistributedAgFEM.jl")

end # module
