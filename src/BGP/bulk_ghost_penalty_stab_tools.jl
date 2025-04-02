"""
    returns the L2 projection in the bounding box space of the trial `u` and test `v` functions after `operation` (e.g. identity or divergence) has been applied. That is, the L2 projection is defined through

        ∫ (operation(u) - Π_{Zbb}(operation(u)) zbb dΩbg_agg_cells = 0   ∀ zbb ∈ Zbb(T), 
        
    with T ∈ T_agg and `Zbb` the bounding box space. Note that Ωbg_agg_cells ⊆ Ωbg_bb. Additionally, Π_{Zbb}(operation(u)) appears on the lhs as the trial function wbb ∈ Zbb.
"""
function setup_L2_proj_in_bb_space(
    dΩbg_agg_cells,             # measure of aggregated cells in background domain
    ref_agg_cell_to_ref_bb_map, # map
    agg_cells_to_aggregate,     # 
    aggregate_to_local_cells,   # 
    u,                          # Trial basis (to project) 
    v,                          # Test basis 
    wbb,                        # Trial basis of bounding box space Zbb
    zbb,                        # Test basis of bounding box space Zbb
    operation,                  # operation to be applied to u and v
    U_agg_cells_local_dof_ids)  # aggregates local dof ids for space U                     
    
    # (0) Obtain triangulation of aggregate cell domain from its measure
    Ωbg_agg_cells=dΩbg_agg_cells.quad.trian

    # (1) Change domain of zbb (test function of Zbb) from Ωbb to Ωbg_agg_cells
    zbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(zbb,
                                                  ref_agg_cell_to_ref_bb_map,
                                                  Ωbg_agg_cells,
                                                  agg_cells_to_aggregate)
    
    # (2) Change domain of wbb (trial function of Zbb) from Ωbb to Ωbg_agg_cells
    wbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(wbb,
                                                  ref_agg_cell_to_ref_bb_map,
                                                  Ωbg_agg_cells,
                                                  agg_cells_to_aggregate)                                               
    # (3) TODO: Compute & assemble contributions to LHS of L2 projection
    agg_cells_to_lhs_contribs=get_array(∫(zbb_Ωbg_agg_cells⋅wbb_Ωbg_agg_cells)dΩbg_agg_cells)
    ass_lhs_map=BulkGhostPenaltyAssembleLhsMap(agg_cells_to_lhs_contribs)
    lhs = lazy_map(ass_lhs_map,aggregate_to_local_cells)

    # (4) Compute & assemble contributions to RHS of L2 projection
    op_u_single_field = operation(_get_single_field_fe_basis(u))
    agg_cells_rhs_contribs=get_array(∫(zbb_Ωbg_agg_cells⋅op_u_single_field)dΩbg_agg_cells)
    ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(U_agg_cells_local_dof_ids,agg_cells_rhs_contribs)
    rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)

    # (5) TO-DO: optimize using our own optimized version Gridap.Fields.Map 
    # of backslash that re-uses storage for lu factors among cells, etc.
    op_v_proj_Zbb_dofs=lazy_map(\,lhs,rhs)

    # (6) Generate bb-wise array of fields. For each aggregate's bounding box,
    # it provides the l2 projection of operation(u) restricted to the cells 
    # included in the bounding box of the aggregate.   
    op_v_proj_Zbb_array=lazy_map(Gridap.Fields.linear_combination,
                                op_v_proj_Zbb_dofs,
                                Gridap.CellData.get_data(zbb))

    # (7) Change domain of proj_op_u and proj_op_v from Ωbb to Ωbg_agg_cells
    op_v_proj_Zbb_array_agg_cells=lazy_map(Broadcasting(∘),
                                           lazy_map(Reindex(op_v_proj_Zbb_array),agg_cells_to_aggregate),
                                           ref_agg_cell_to_ref_bb_map)
    op_u_proj_Zbb_array_agg_cells=lazy_map(transpose, op_v_proj_Zbb_array_agg_cells)
    if (_is_multifield_fe_basis_component(u))
        @assert _is_multifield_fe_basis_component(v)
        @assert _nfields(u)==_nfields(v)
        nfields=_nfields(u)
        fieldid=_fieldid(u)
        op_u_proj_Zbb_array_agg_cells=lazy_map(
                                    Gridap.Fields.BlockMap((1,nfields),fieldid),
                                    op_u_proj_Zbb_array_agg_cells)
    end
    op_u_proj_Zbb_agg_cells = Gridap.CellData.GenericCellField(op_u_proj_Zbb_array_agg_cells,Ωbg_agg_cells,ReferenceDomain()) 

    if (_is_multifield_fe_basis_component(v))
            @assert _is_multifield_fe_basis_component(v)
            @assert _nfields(u)==_nfields(v)
            nfields=_nfields(v)
            fieldid=_fieldid(v)
            op_v_proj_Zbb_array_agg_cells=lazy_map(
                                        Gridap.Fields.BlockMap(nfields,fieldid),
                                        op_v_proj_Zbb_array_agg_cells)
    end 

    op_v_proj_Zbb_agg_cells = Gridap.CellData.GenericCellField(op_v_proj_Zbb_array_agg_cells,Ωbg_agg_cells,ReferenceDomain()) 

    op_u_proj_Zbb_agg_cells, op_v_proj_Zbb_agg_cells
end

# TODO: the following is not working yet:
# """
#     returns the L2 projection in the bounding box space of a function `f`. That is, the L2 projection is defined through

#         ∫ (f - Π_{Zbb}(f) zbb dΩbg_agg_cells = 0   ∀ zbb ∈ Zbb(T), 
        
#     with T ∈ T_agg and `Zbb` the bounding box space. Note that Ωbg_agg_cells ⊆ Ωbg_bb. Additionally, Π_{Zbb}(f) appears on the lhs as the trial function wbb ∈ Zbb.
# """
# function setup_L2_proj_in_bb_space(
#     dΩbg_agg_cells,             # measure of aggregated cells in background domain
#     ref_agg_cell_to_ref_bb_map, # map
#     agg_cells_to_aggregate,     # 
#     aggregate_to_local_cells,   # 
#     f,
#     u,                          # Trial basis (to project) 
#     v,                          # Test basis 
#     wbb,                        # Trial basis of bounding box space Zbb
#     zbb,                        # Test basis of bounding box space Zbb
#     operation,                  # operation to be applied to u and v
#     U_agg_cells_local_dof_ids)  # aggregates local dof ids for space U                     
    
#     # (0) Obtain triangulation of aggregate cell domain from its measure
#     Ωbg_agg_cells=dΩbg_agg_cells.quad.trian

#     # (1) Change domain of zbb (test function of Zbb) from Ωbb to Ωbg_agg_cells
#     zbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(zbb,
#                                                   ref_agg_cell_to_ref_bb_map,
#                                                   Ωbg_agg_cells,
#                                                   agg_cells_to_aggregate)
    
#     # (2) Change domain of wbb (trial function of Zbb) from Ωbb to Ωbg_agg_cells
#     wbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(wbb,
#                                                   ref_agg_cell_to_ref_bb_map,
#                                                   Ωbg_agg_cells,
#                                                   agg_cells_to_aggregate)                                               
#     # (3) TODO: Compute & assemble contributions to LHS of L2 projection
#     agg_cells_to_lhs_contribs=get_array(∫(zbb_Ωbg_agg_cells⋅wbb_Ωbg_agg_cells)dΩbg_agg_cells)
#     ass_lhs_map=BulkGhostPenaltyAssembleLhsMap(agg_cells_to_lhs_contribs)
#     lhs = lazy_map(ass_lhs_map,aggregate_to_local_cells)

#     # (4) Compute & assemble contributions to RHS of L2 projection
#     op_u_single_field = operation(_get_single_field_fe_basis(u))
#     agg_cells_rhs_contribs=get_array(∫(zbb_Ωbg_agg_cells⋅op_u_single_field)dΩbg_agg_cells)
#     ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(U_agg_cells_local_dof_ids,agg_cells_rhs_contribs)
#     rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)

#     # (5) TO-DO: optimize using our own optimized version Gridap.Fields.Map 
#     # of backslash that re-uses storage for lu factors among cells, etc.
#     op_v_proj_Zbb_dofs=lazy_map(\,lhs,rhs)

#     # (6) Generate bb-wise array of fields. For each aggregate's bounding box,
#     # it provides the l2 projection of operation(u) restricted to the cells 
#     # included in the bounding box of the aggregate.   
#     op_v_proj_Zbb_array=lazy_map(Gridap.Fields.linear_combination,
#                                 op_v_proj_Zbb_dofs,
#                                 Gridap.CellData.get_data(zbb))

#     # (7) Change domain of proj_op_u and proj_op_v from Ωbb to Ωbg_agg_cells
#     op_v_proj_Zbb_array_agg_cells=lazy_map(Broadcasting(∘),
#                                            lazy_map(Reindex(op_v_proj_Zbb_array),agg_cells_to_aggregate),
#                                            ref_agg_cell_to_ref_bb_map)
#     op_u_proj_Zbb_array_agg_cells=lazy_map(transpose, op_v_proj_Zbb_array_agg_cells)
#     if (_is_multifield_fe_basis_component(u))
#         @assert _is_multifield_fe_basis_component(v)
#         @assert _nfields(u)==_nfields(v)
#         nfields=_nfields(u)
#         fieldid=_fieldid(u)
#         op_u_proj_Zbb_array_agg_cells=lazy_map(
#                                     Gridap.Fields.BlockMap((1,nfields),fieldid),
#                                     op_u_proj_Zbb_array_agg_cells)
#     end
#     op_u_proj_Zbb_agg_cells = Gridap.CellData.GenericCellField(op_u_proj_Zbb_array_agg_cells,Ωbg_agg_cells,ReferenceDomain()) 

#     if (_is_multifield_fe_basis_component(v))
#             @assert _is_multifield_fe_basis_component(v)
#             @assert _nfields(u)==_nfields(v)
#             nfields=_nfields(v)
#             fieldid=_fieldid(v)
#             op_v_proj_Zbb_array_agg_cells=lazy_map(
#                                         Gridap.Fields.BlockMap(nfields,fieldid),
#                                         op_v_proj_Zbb_array_agg_cells)
#     end 

#     op_v_proj_Zbb_agg_cells = Gridap.CellData.GenericCellField(op_v_proj_Zbb_array_agg_cells,Ωbg_agg_cells,ReferenceDomain()) 

#     op_u_proj_Zbb_agg_cells, op_v_proj_Zbb_agg_cells
# end

"""
  Compute and assemble the bulk penalty stabilization term

  γ ∫( (operation_u(u) - proj_op_u)⊙(operation_v(v) - proj_op_v) )*dD

  with `u` and `v` the trial and test basis functions (which may be components of a MultiField) and do not have to belong to the same space. The operations `operation_u` and `operation_v` (e.g. identity or divergence) are applied to u and v, respectively. The L2 projections `proj_op_u` and `proj_op_v` can be computed via `setup_L2_proj_in_bb_space` and are to be passed as arguments. `dD` is the measure of the cut cells or full aggregate, that is `dΩbg_cut_cells` or `dΩbg_agg_cells`.
"""
function bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                                                                    Ωbg_agg_cells,
                                                                    γ,  
                                                                    u, 
                                                                    v, 
                                                                    proj_op_u,
                                                                    proj_op_v,
                                                                    dof_ids_u,
                                                                    dof_ids_v,
                                                                    dof_ids_proj_u,
                                                                    dof_ids_proj_v,
                                                                    operation_u,
                                                                    operation_v)

    # # BlockMap preparatory steps for the dof ids
    # if (_is_multifield_fe_basis_component(u))
    #     nfields=_nfields(u)
    #     fieldid=_fieldid(u)
    #     dof_ids_u=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_u)
    # end
    # if (_is_multifield_fe_basis_component(v))
    #     nfields=_nfields(v)
    #     fieldid=_fieldid(v)
    #     dof_ids_v=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_v)
    # end
    # if (_is_multifield_fe_basis_component(u))
    #     nfields=_nfields(u)
    #     fieldid=_fieldid(u)
    #     dof_ids_proj_u=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_proj_u)
    # end
    # if (_is_multifield_fe_basis_component(v))
    #     nfields=_nfields(v)
    #     fieldid=_fieldid(v)
    #     dof_ids_proj_v=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_proj_v)
    # end

    # Manually set up the arrays that collect_cell_matrix would return automatically
    w = []
    r = []
    c = []

    # (1) op_u⊙op_v term
    op_u_op_v_mat_contribs=get_array(∫(γ*operation_u(u)⊙operation_v(v))*dD)
    push!(w, op_u_op_v_mat_contribs)
    push!(r, dof_ids_v) 
    push!(c, dof_ids_u)

    # (2) proj_op_u⊙op_v term

    # In the MultiField case, I have had to add this change domain 
    # call before setting up the term right below. Otherwise, we get an error 
    # when trying to multiply the fields. Not sure why this is happening
    op_v_on_Ωbg_agg_cells = Gridap.CellData.change_domain(operation_v(v),Ωbg_agg_cells,ReferenceDomain())

    proj_op_u_op_v_mat_contribs=get_array(∫(γ*(-1.0)*(proj_op_u⊙op_v_on_Ωbg_agg_cells))*dD)
    push!(w, proj_op_u_op_v_mat_contribs)    
    push!(r, dof_ids_v) 
    push!(c, dof_ids_proj_u)

    # (3) proj_op_u⊙proj_op_v
    proj_op_u_proj_op_v_mat_contribs=get_array(∫(γ*(proj_op_u⊙proj_op_v))*dD)
    push!(w, proj_op_u_proj_op_v_mat_contribs)
    push!(r, dof_ids_proj_v) 
    push!(c, dof_ids_proj_u)

    # (4) op_u⊙proj_op_v

    # In the MultiField case, I have had to add this change domain 
    # call before setting up the term right below. Otherwise, we get an error 
    # when trying to multiply the fields. Not sure why this is happening
    op_u_on_Ωbg_agg_cells = Gridap.CellData.change_domain(operation_u(u),Ωbg_agg_cells,ReferenceDomain())

    op_u_proj_op_v_mat_contribs=get_array(∫(-γ*op_u_on_Ωbg_agg_cells⊙proj_op_v)*dD)
    push!(w, op_u_proj_op_v_mat_contribs)
    push!(r, dof_ids_proj_v) 
    push!(c, dof_ids_u)

    w, r, c
end

"""
  Compute and assemble the bulk penalty stabilization term

  γ ∫( (operation_u(u) - proj_op_u)⊙(operation_v(v)) )*dD

  with `u` and `v` the trial and test basis functions (which may be components of a MultiField). The operations `operation_u` and `operation_v` (e.g. identity or divergence) are applied to u and v, respectively. The L2 projections `proj_op_u` and `proj_op_v` can be computed via `setup_L2_proj_in_bb_space` and are to be passed as arguments. `dD` is the measure of the cut cells or full aggregate, that is `dΩbg_cut_cells` or `dΩbg_agg_cells`.

"""
function bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dD,
                                                                    Ωbg_agg_cells,
                                                                    γ,  
                                                                    u, 
                                                                    v, 
                                                                    proj_op_u,
                                                                    dof_ids_u,
                                                                    dof_ids_v,
                                                                    dof_ids_proj_u,
                                                                    operation_u,
                                                                    operation_v)

    # # BlockMap preparatory steps for the dof ids
    # if (_is_multifield_fe_basis_component(u))
    #     nfields=_nfields(u)
    #     fieldid=_fieldid(u)
    #     dof_ids_u=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_u)
    # end
    # if (_is_multifield_fe_basis_component(v))
    #     nfields=_nfields(v)
    #     fieldid=_fieldid(v)
    #     dof_ids_v=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_v)
    # end
    # if (_is_multifield_fe_basis_component(u))
    #     nfields=_nfields(u)
    #     fieldid=_fieldid(u)
    #     dof_ids_proj_u=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_proj_u)
    # end

    # Manually set up the arrays that collect_cell_matrix would return automatically
    w = []
    r = []
    c = []

    # (1) op_u⊙op_v term
    op_u_op_v_mat_contribs=get_array(∫(γ*operation_u(u)⊙operation_v(v))*dD)
    push!(w, op_u_op_v_mat_contribs)
    push!(r, dof_ids_v) 
    push!(c, dof_ids_u)

    # (2) proj_op_u⊙op_v term

    # In the MultiField case, I have had to add this change domain 
    # call before setting up the term right below. Otherwise, we get an error 
    # when trying to multiply the fields. Not sure why this is happening
    op_v_on_Ωbg_agg_cells = Gridap.CellData.change_domain(operation_v(v),Ωbg_agg_cells,ReferenceDomain())

    proj_op_u_op_v_mat_contribs=get_array(∫(γ*(-1.0)*(proj_op_u⊙op_v_on_Ωbg_agg_cells))*dD)
    push!(w, proj_op_u_op_v_mat_contribs)    
    push!(r, dof_ids_v) 
    push!(c, dof_ids_proj_u)

    w, r, c
end

"""
  Compute and assemble the bulk penalty stabilization term

  γ ∫( (operation_v(v) - proj_op_v)⊙(rhs_function) )*dD

  with `v` the test basis functions (which may be components of a MultiField). The operation `operation_v` (e.g. identity or divergence) is applied to v. The L2 projections `proj_op_u` can be computed via `setup_L2_proj_in_bb_space`. `dD` is the measure of the cut cells or full aggregate, that is `dΩbg_cut_cells` or `dΩbg_agg_cells`. The function `rhs_function` is the forcing term appearing on the right hand-side of the PDE. 

"""
function bulk_ghost_penalty_stabilization_collect_cell_vector_on_D(dD,
                                                                    γ,  
                                                                    v, 
                                                                    proj_op_v,
                                                                    dof_ids_v,
                                                                    dof_ids_proj_v,
                                                                    operation_v,
                                                                    rhs_function)
    # # BlockMap preparatory steps for the dof ids
    # if (_is_multifield_fe_basis_component(v))
    #     nfields=_nfields(v)
    #     fieldid=_fieldid(v)
    #     dof_ids_v=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_v)
    # end
    # if (_is_multifield_fe_basis_component(v))
    #     nfields=_nfields(v)
    #     fieldid=_fieldid(v)
    #     dof_ids_proj_v=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_proj_v)
    # end

    # Manually set up the arrays that collect_cell_vector would return automatically
    w = []
    r = []

    # (1) op_v⊙rhs_function term
    op_v_rhs_vec_contribs=get_array(∫(γ*operation_v(v)⊙(rhs_function))*dD)
    push!(w, op_v_rhs_vec_contribs)
    push!(r, dof_ids_v) 

    # (2) proj_op_v⊙rhs_function
    proj_op_v_rhs_vec_contribs=get_array(∫(γ*(-1.0)*(proj_op_v⊙(rhs_function)))*dD)
    push!(w, proj_op_v_rhs_vec_contribs)    
    push!(r, dof_ids_proj_v) 

    w, r
end

"""
  Compute and assemble the bulk penalty stabilization term

  γ ∫( (operation_v(v))⊙(rhs_function) )*dD

  with `v` the test basis functions (which may be components of a MultiField). The operation `operation_v` (e.g. identity or divergence) is applied to v. The L2 projections `proj_op_u` can be computed via `setup_L2_proj_in_bb_space`. `dD` is the measure of the cut cells or full aggregate, that is `dΩbg_cut_cells` or `dΩbg_agg_cells`. The function `rhs_function` is the forcing term appearing on the right hand-side of the PDE. 

"""
function bulk_ghost_penalty_stabilization_collect_cell_vector_on_D(dD,
                                                                    γ,  
                                                                    v, 
                                                                    dof_ids_v,
                                                                    operation_v,
                                                                    rhs_function)
    # # BlockMap preparatory steps for the dof ids
    # if (_is_multifield_fe_basis_component(v))
    #     nfields=_nfields(v)
    #     fieldid=_fieldid(v)
    #     dof_ids_v=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_v)
    # end
    # if (_is_multifield_fe_basis_component(v))
    #     nfields=_nfields(v)
    #     fieldid=_fieldid(v)
    #     dof_ids_proj_v=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),dof_ids_proj_v)
    # end

    # Manually set up the arrays that collect_cell_vector would return automatically
    w = []
    r = []

    # (1) op_v⊙rhs_function term
    op_v_rhs_vec_contribs=get_array(∫(γ*operation_v(v)⊙(rhs_function))*dD)
    push!(w, op_v_rhs_vec_contribs)
    push!(r, dof_ids_v) 

    w, r
end