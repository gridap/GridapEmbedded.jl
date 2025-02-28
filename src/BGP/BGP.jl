using Gridap
using FillArrays
using LinearAlgebra

include("aggregates_bounding_boxes_tools.jl")

export setup_aggregate_to_cells
export setup_aggregates_bounding_box_model
export flatten
export restrict_cells
export setup_cells_to_aggregate
export setup_ref_agg_cell_to_ref_bb_map
export setup_aggregate_to_local_cells
export compute_agg_cells_local_dof_ids
export compute_aggregate_dof_ids

include("bulk_ghost_penalty_stab_tools.jl")
export set_up_bulk_ghost_penalty_lhs
export bulk_ghost_penalty_stabilization_collect_cell_matrix_on_cut_cells, bulk_ghost_penalty_stabilization_collect_cell_matrix_on_cut_cells_full
export div_penalty_stabilization_collect_cell_matrix_on_cut_cells, div_penalty_stabilization_collect_cell_matrix_on_cut_cells_full
export pmix_penalty_stabilization_collect_cell_matrix_on_cut_cells, dmix_penalty_stabilization_collect_cell_matrix_on_cut_cells, dmix_penalty_stabilization_collect_cell_vector_on_cut_cells

include("fields_and_blocks_tools.jl")
export _restrict_to_block
# export _restrict_to_block, _get_single_field_fe_basis, _get_single_field_fe_basis,  _is_multifield_fe_basis_component,  _is_multifield_fe_basis_component, _nfields, _fieldid

include("BulkGhostPenaltyAssembleMaps.jl")

# include("darcy_preconditioning_tools.jl")
# export assemble_darcy_preconditioner_matrix


