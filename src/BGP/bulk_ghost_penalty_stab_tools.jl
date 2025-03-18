"""
  dv, du: Test and trial basis functions. They may be components of a MultiFieldCellField

  # Compute and assemble the bulk penalty stabilization term
  # ∫( (dv-dv_l2_proj_agg_cells)*(du-du_l2_proj_agg_cells))*dΩ_agg_cells (long version)
  # ∫( (dv)*(du-du_l2_proj_agg_cells))*dΩ_agg_cells (simplified, equivalent version)
"""
function bulk_ghost_penalty_stabilization_collect_cell_matrix(aggregate_to_local_cells,
                                                              agg_cells_to_aggregate,
                                                              ref_agg_cell_to_ref_bb_map,
                                                              dΩagg_cells,
                                                              dv,    # Test basis
                                                              du,    # Trial basis (to project)
                                                              dvbb,  # Bounding box space test basis
                                                              lhs,
                                                              Ωagg_cell_dof_ids,
                                                              agg_cells_local_dof_ids,
                                                              agg_cells_to_aggregate_dof_ids,
                                                              γ)


    Ωagg_cells=dΩagg_cells.quad.trian

    # Change domain of vbb (test) from Ωbb to Ωagg_cells
    dvbb_Ωagg_cells=change_domain_bb_to_agg_cells(dvbb,
                                                  ref_agg_cell_to_ref_bb_map,
                                                  Ωagg_cells,
                                                  agg_cells_to_aggregate)             

    du_single_field=_get_single_field_fe_basis(du)                                              
    agg_cells_rhs_contribs=get_array(∫(dvbb_Ωagg_cells⋅du_single_field)dΩagg_cells)
    ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(agg_cells_local_dof_ids,agg_cells_rhs_contribs)
    rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)

    # TO-DO: optimize using our own optimized version Gridap.Fields.Map 
    # of backslash that re-uses storage for lu factors among cells, etc.
    dv_l2_proj_bb_dofs=lazy_map(\,lhs,rhs)

    # Generate bb-wise array of fields. For each aggregate's bounding box,
    # it provides the l2 projection of all basis functions in du 
    # restricted to the cells included in the bounding box of the aggregate   
    dv_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination,
                                dv_l2_proj_bb_dofs,
                                Gridap.CellData.get_data(dvbb))

    # # Change domain of dv_l2_proj_bb_array from bb to agg_cells
    dv_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),
                                           lazy_map(Reindex(dv_l2_proj_bb_array),agg_cells_to_aggregate),
                                           ref_agg_cell_to_ref_bb_map)

    du_l2_proj_bb_array_agg_cells=lazy_map(transpose, dv_l2_proj_bb_array_agg_cells)
    
    if (_is_multifield_fe_basis_component(du))
        @assert _is_multifield_fe_basis_component(dv)
        @assert _nfields(du)==_nfields(dv)
        nfields=_nfields(du)
        fieldid=_fieldid(du)
        du_l2_proj_bb_array_agg_cells=lazy_map(
                                    Gridap.Fields.BlockMap((1,nfields),fieldid),
                                    du_l2_proj_bb_array_agg_cells)
    end 

    du_l2_proj_agg_cells = Gridap.CellData.GenericCellField(du_l2_proj_bb_array_agg_cells,
                                                            dΩagg_cells.quad.trian,
                                                            ReferenceDomain())

    # Manually set up the arrays that collect_cell_matrix would return automatically
    w = []
    r = []
    c = []

    dv_du_mat_contribs=get_array(∫(γ*dv⋅du)*dΩagg_cells)
    if (_is_multifield_fe_basis_component(dv))
        nfields=_nfields(dv)
        fieldid=_fieldid(dv)
        Ωagg_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),Ωagg_cell_dof_ids)
    end
    push!(w, dv_du_mat_contribs)
    push!(r, Ωagg_cell_dof_ids)
    push!(c, Ωagg_cell_dof_ids)

    # In the MultiField case, I have had to add this change domain 
    # call before setting up the term right below. Otherwise, we get an error 
    # when trying to multiply the fields. Not sure why this is happening
    dvΩagg_cells=Gridap.CellData.change_domain(dv,Ωagg_cells,ReferenceDomain())

    proj_dv_du_mat_contribs=get_array(∫(γ*(-1.0)*dvΩagg_cells⋅(du_l2_proj_agg_cells))*dΩagg_cells)
    
    if (_is_multifield_fe_basis_component(du))
        nfields=_nfields(du)
        fieldid=_fieldid(du)
        agg_cells_to_aggregate_dof_ids=
           lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),agg_cells_to_aggregate_dof_ids)
    end
    
    push!(w, proj_dv_du_mat_contribs)
    push!(r, Ωagg_cell_dof_ids)
    push!(c, agg_cells_to_aggregate_dof_ids)

    w, r, c
end

## Create collect function on cut cells only
"""
  dv, du: Test and trial basis functions. They may be components of a MultiFieldCellField

  # Compute and assemble the bulk penalty stabilization term
  # ∫( (dv)*(du-du_l2_proj_agg_cells))*dΩ_cut_cells
"""
function bulk_ghost_penalty_stabilization_collect_cell_matrix_on_cut_cells(aggregate_to_local_cells,
                                                              agg_cells_to_aggregate,
                                                              ref_agg_cell_to_ref_bb_map,
                                                              dΩbg_agg_cells,
                                                              dv,    # Test basis
                                                              du,    # Trial basis (to project)
                                                              dvbb,  # Bounding box space test basis
                                                              lhs,
                                                              Ωbg_cut_cell_dof_ids,
                                                              agg_cells_local_dof_ids,
                                                              cut_cells_to_aggregate_dof_ids,
                                                              γ,
                                                              dΩbg_cut_cells) 


    Ωbg_agg_cells=dΩbg_agg_cells.quad.trian

    # Change domain of vbb (test) from Ωbb to Ωagg_cells
    dvbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(dvbb,
                                                  ref_agg_cell_to_ref_bb_map,
                                                  Ωbg_agg_cells,
                                                  agg_cells_to_aggregate)             

    du_single_field=_get_single_field_fe_basis(du)                                              
    agg_cells_rhs_contribs=get_array(∫(dvbb_Ωbg_agg_cells⋅du_single_field)dΩbg_agg_cells)
    ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(agg_cells_local_dof_ids,agg_cells_rhs_contribs)
    rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)

    # TO-DO: optimize using our own optimized version Gridap.Fields.Map 
    # of backslash that re-uses storage for lu factors among cells, etc.
    dv_l2_proj_bb_dofs=lazy_map(\,lhs,rhs)

    # Generate bb-wise array of fields. For each aggregate's bounding box,
    # it provides the l2 projection of all basis functions in du 
    # restricted to the cells included in the bounding box of the aggregate   
    dv_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination,
                                dv_l2_proj_bb_dofs,
                                Gridap.CellData.get_data(dvbb))

    # # Change domain of dv_l2_proj_bb_array from bb to agg_cells
    dv_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),
                                           lazy_map(Reindex(dv_l2_proj_bb_array),agg_cells_to_aggregate),
                                           ref_agg_cell_to_ref_bb_map)

    du_l2_proj_bb_array_agg_cells=lazy_map(transpose, dv_l2_proj_bb_array_agg_cells)
    
    if (_is_multifield_fe_basis_component(du))
        @assert _is_multifield_fe_basis_component(dv)
        @assert _nfields(du)==_nfields(dv)
        nfields=_nfields(du)
        fieldid=_fieldid(du)
        du_l2_proj_bb_array_agg_cells=lazy_map(
                                    Gridap.Fields.BlockMap((1,nfields),fieldid),
                                    du_l2_proj_bb_array_agg_cells)
    end 

    du_l2_proj_agg_cells = Gridap.CellData.GenericCellField(du_l2_proj_bb_array_agg_cells, Ωbg_agg_cells,ReferenceDomain()) 

    # Manually set up the arrays that collect_cell_matrix would return automatically
    w = []
    r = []
    c = []

    dv_du_mat_contribs=get_array(∫(γ*dv⋅du)*dΩbg_cut_cells)

    if (_is_multifield_fe_basis_component(dv))
        nfields=_nfields(dv)
        fieldid=_fieldid(dv)
        Ωbg_cut_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),Ωbg_cut_cell_dof_ids)
    end
    push!(w, dv_du_mat_contribs)
    push!(r, Ωbg_cut_cell_dof_ids) 
    push!(c, Ωbg_cut_cell_dof_ids)

    # In the MultiField case, I have had to add this change domain 
    # call before setting up the term right below. Otherwise, we get an error 
    # when trying to multiply the fields. Not sure why this is happening
    dvΩagg_cells=Gridap.CellData.change_domain(dv,Ωbg_agg_cells,ReferenceDomain())

    proj_dv_du_mat_contribs=get_array(∫(γ*(-1.0)*dvΩagg_cells⋅(du_l2_proj_agg_cells))*dΩbg_cut_cells)

    if (_is_multifield_fe_basis_component(du))
        nfields=_nfields(du)
        fieldid=_fieldid(du)
        cut_cells_to_aggregate_dof_ids=
        lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),cut_cells_to_aggregate_dof_ids)
    end

    push!(w, proj_dv_du_mat_contribs)    
    push!(r, Ωbg_cut_cell_dof_ids) 
    push!(c, cut_cells_to_aggregate_dof_ids)

    w, r, c
end

## Create collect function on cut cells only
"""
  dv, du: Test and trial basis functions. They may be components of a MultiFieldCellField

  # Compute and assemble the bulk penalty stabilization term
  # ∫( (dv-dv_l2_proj_agg_cells)*(du-du_l2_proj_agg_cells))*dΩ_cut_cells
"""
function bulk_ghost_penalty_stabilization_collect_cell_matrix_on_cut_cells_full(aggregate_to_local_cells,
                                                              agg_cells_to_aggregate,
                                                              ref_agg_cell_to_ref_bb_map,
                                                              dΩbg_agg_cells,
                                                              dv,    # Test basis
                                                              du,    # Trial basis (to project)
                                                              dvbb,  # Bounding box space test basis
                                                              lhs,
                                                              Ωbg_cut_cell_dof_ids,
                                                              agg_cells_local_dof_ids,
                                                              cut_cells_to_aggregate_dof_ids,
                                                              γ,
                                                              dΩbg_cut_cells) 


    Ωbg_agg_cells=dΩbg_agg_cells.quad.trian

    ## (I) Compute projections dv_l2_proj_agg_cells & du_l2_proj_agg_cells
    # Change domain of vbb (test) from Ωbb to Ωagg_cells
    dvbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(dvbb,
                                                  ref_agg_cell_to_ref_bb_map,
                                                  Ωbg_agg_cells,
                                                  agg_cells_to_aggregate)             

    du_single_field=_get_single_field_fe_basis(du)                                              
    agg_cells_rhs_contribs=get_array(∫(dvbb_Ωbg_agg_cells⋅du_single_field)dΩbg_agg_cells)
    ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(agg_cells_local_dof_ids,agg_cells_rhs_contribs)
    rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)

    # TO-DO: optimize using our own optimized version Gridap.Fields.Map 
    # of backslash that re-uses storage for lu factors among cells, etc.
    dv_l2_proj_bb_dofs=lazy_map(\,lhs,rhs)

    # Generate bb-wise array of fields. For each aggregate's bounding box,
    # it provides the l2 projection of all basis functions in du 
    # restricted to the cells included in the bounding box of the aggregate   
    dv_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination,
                                dv_l2_proj_bb_dofs,
                                Gridap.CellData.get_data(dvbb))

    # # Change domain of dv_l2_proj_bb_array from bb to agg_cells
    dv_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),
                                           lazy_map(Reindex(dv_l2_proj_bb_array),agg_cells_to_aggregate),
                                           ref_agg_cell_to_ref_bb_map)

    du_l2_proj_bb_array_agg_cells=lazy_map(transpose, dv_l2_proj_bb_array_agg_cells)
    
    if (_is_multifield_fe_basis_component(du))
        @assert _is_multifield_fe_basis_component(dv)
        @assert _nfields(du)==_nfields(dv)
        nfields=_nfields(du)
        fieldid=_fieldid(du)
        du_l2_proj_bb_array_agg_cells=lazy_map(
                                    Gridap.Fields.BlockMap((1,nfields),fieldid),
                                    du_l2_proj_bb_array_agg_cells)
    end

    if (_is_multifield_fe_basis_component(dv))
        @assert _is_multifield_fe_basis_component(dv)
        @assert _nfields(du)==_nfields(dv)
        nfields=_nfields(dv)
        fieldid=_fieldid(dv)
        dv_l2_proj_bb_array_agg_cells=lazy_map(
                                    Gridap.Fields.BlockMap(nfields,fieldid),
                                    dv_l2_proj_bb_array_agg_cells)
     end 

    du_l2_proj_agg_cells = Gridap.CellData.GenericCellField(du_l2_proj_bb_array_agg_cells, Ωbg_agg_cells,ReferenceDomain()) 
    dv_l2_proj_agg_cells = Gridap.CellData.GenericCellField(dv_l2_proj_bb_array_agg_cells, Ωbg_agg_cells,ReferenceDomain()) 

    # BlockMap preparatory steps for the dof ids
    if (_is_multifield_fe_basis_component(dv))
        nfields=_nfields(dv)
        fieldid=_fieldid(dv)
        Ωbg_cut_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),Ωbg_cut_cell_dof_ids)
    end

    if (_is_multifield_fe_basis_component(du))
        nfields=_nfields(du)
        fieldid=_fieldid(du)
        cut_cells_to_aggregate_dof_ids=
        lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),cut_cells_to_aggregate_dof_ids)
    end

    # Manually set up the arrays that collect_cell_matrix would return automatically
    w = []
    r = []
    c = []

    ## (1) dv*du term
    dv_du_mat_contribs=get_array(∫(γ*dv⋅du)*dΩbg_cut_cells)

    push!(w, dv_du_mat_contribs)
    push!(r, Ωbg_cut_cell_dof_ids) 
    push!(c, Ωbg_cut_cell_dof_ids)

    ## (2) dv*proj_du term
    # In the MultiField case, I have had to add this change domain 
    # call before setting up the term right below. Otherwise, we get an error 
    # when trying to multiply the fields. Not sure why this is happening
    dvΩagg_cells=Gridap.CellData.change_domain(dv,Ωbg_agg_cells,ReferenceDomain())

    dv_proj_du_mat_contribs=get_array(∫(γ*(-1.0)*dvΩagg_cells⋅(du_l2_proj_agg_cells))*dΩbg_cut_cells)

    push!(w, dv_proj_du_mat_contribs)    
    push!(r, Ωbg_cut_cell_dof_ids) 
    push!(c, cut_cells_to_aggregate_dof_ids)

    ## (3) proj_dv*proj_du
    proj_dv_proj_du_mat_contribs=get_array(∫(γ*dv_l2_proj_agg_cells⋅(du_l2_proj_agg_cells))*dΩbg_cut_cells)

    push!(w, proj_dv_proj_du_mat_contribs)
    push!(r, cut_cells_to_aggregate_dof_ids) 
    push!(c, cut_cells_to_aggregate_dof_ids)

    ## (4) proj_dv*du
    duΩagg_cells=Gridap.CellData.change_domain(du,Ωbg_agg_cells,ReferenceDomain())
    proj_dv_du_mat_contribs=get_array(∫(-γ*dv_l2_proj_agg_cells⋅duΩagg_cells)*dΩbg_cut_cells)

    push!(w, proj_dv_du_mat_contribs)
    push!(r, cut_cells_to_aggregate_dof_ids) 
    push!(c, Ωbg_cut_cell_dof_ids)

    w, r, c
end

function div_penalty_stabilization_collect_cell_matrix(aggregate_to_local_cells,
                                                       agg_cells_to_aggregate,
                                                       ref_agg_cell_to_ref_bb_map,
                                                       dΩagg_cells,
                                                       dv,    # Test basis
                                                       du,    # Trial basis (to project)
                                                       qbb,   # Bounding box space test basis
                                                       p_lhs,
                                                       U_Ωagg_cell_dof_ids, 
                                                       U_agg_cells_local_dof_ids,
                                                       U_agg_cells_to_aggregate_dof_ids,
                                                       γ)
    Ωagg_cells=dΩagg_cells.quad.trian

    ## +DIVU STABILIZATION 
    ## Compute ∫( (div_dv)*(div_du-div_du_l2_proj_agg_cells))*dΩ_agg_cells (simplified, equivalent version)
    # Manually set up the arrays that collect_cell_matrix would return automatically
    w = []
    r = []
    c = []

    if (_is_multifield_fe_basis_component(dv))
        nfields=_nfields(dv)
        fieldid=_fieldid(dv)
        U_Ωagg_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),U_Ωagg_cell_dof_ids)
    end

    ## (I) Let's set up the div_dv_div_du_mat_contribs 
    # (first part of the integral, that is ∫( (div_dv)*(div_du)*dΩ_agg_cells
    div_dv_div_du_mat_contribs=get_array(∫(γ*(∇⋅dv)⋅(∇⋅du))*dΩagg_cells)
    push!(w, div_dv_div_du_mat_contribs)
    push!(r, U_Ωagg_cell_dof_ids)
    push!(c, U_Ωagg_cell_dof_ids)
    div_du=∇⋅(_get_single_field_fe_basis(du))

    ### Compute Π_Q_bb(div_du)
    dqbb_Ωagg_cells =change_domain_bb_to_agg_cells(qbb,ref_agg_cell_to_ref_bb_map,Ωagg_cells,agg_cells_to_aggregate)  
    agg_cells_rhs_contribs=get_array(∫(dqbb_Ωagg_cells⋅div_du)dΩagg_cells) 
    ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(U_agg_cells_local_dof_ids,agg_cells_rhs_contribs)
    rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)
    div_dv_l2_proj_bb_dofs=lazy_map(\,p_lhs,rhs)

    # Generate bb-wise array of fields. For each aggregate's bounding box,
        # it provides the l2 projection of all basis functions in dp 
        # restricted to the cells included in the bounding box of the aggregate  
    div_dv_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination,
                                    div_dv_l2_proj_bb_dofs,
                                    Gridap.CellData.get_data(qbb))

    div_dv_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),
        lazy_map(Reindex(div_dv_l2_proj_bb_array),agg_cells_to_aggregate),
        ref_agg_cell_to_ref_bb_map)

    div_du_l2_proj_bb_array_agg_cells=lazy_map(transpose, div_dv_l2_proj_bb_array_agg_cells)

    if (_is_multifield_fe_basis_component(du))
        @assert _is_multifield_fe_basis_component(du)
        @assert _nfields(du)==_nfields(dv)
        nfields=_nfields(du)
        fieldid=_fieldid(dv)
        div_du_l2_proj_bb_array_agg_cells=lazy_map(
                                    Gridap.Fields.BlockMap((fieldid,nfields),fieldid),
                                    div_du_l2_proj_bb_array_agg_cells)
    end 

    div_du_l2_proj_agg_cells = 
        Gridap.CellData.GenericCellField(div_du_l2_proj_bb_array_agg_cells, dΩagg_cells.quad.trian,ReferenceDomain())

    # In the MultiField case, I have had to add this change domain 
    # call before setting up the term right below. Otherwise, we get an error 
    # when trying to multiply the fields. Not sure why this is happening
    # dv_Ωagg_cells=Gridap.CellData.change_domain(dv,Ωagg_cells,ReferenceDomain())
    div_dv_Ωagg_cells=Gridap.CellData.change_domain(∇⋅dv,Ωagg_cells,ReferenceDomain()) 

    proj_div_dv_div_du_mat_contribs=get_array(∫(γ*(-1.0)*div_dv_Ωagg_cells*(div_du_l2_proj_agg_cells))*dΩagg_cells)

    if (_is_multifield_fe_basis_component(du))
        nfields=_nfields(du)
        fieldid=_fieldid(du)
        U_agg_cells_to_aggregate_dof_ids=
           lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),U_agg_cells_to_aggregate_dof_ids)
    end

    push!(w, proj_div_dv_div_du_mat_contribs)
    push!(r, U_Ωagg_cell_dof_ids)              # test
    push!(c, U_agg_cells_to_aggregate_dof_ids) # trial 
    w,r,c
end


function div_penalty_stabilization_collect_cell_matrix_on_cut_cells(aggregate_to_local_cells,
    agg_cells_to_aggregate,
    ref_agg_cell_to_ref_bb_map,
    dΩbg_agg_cells,
    dv,    # Test basis
    du,    # Trial basis (to project)
    qbb,   # Bounding box space test basis
    p_lhs,
    U_Ωbg_cut_cell_dof_ids, 
    U_agg_cells_local_dof_ids,
    U_cut_cells_to_aggregate_dof_ids,
    γ,
    dΩbg_cut_cells)
    
    Ωbg_agg_cells=dΩbg_agg_cells.quad.trian
 
    ## Compute ∫( (div_dv)*(div_du-div_du_l2_proj_agg_cells))*dΩbg_cut_cells 
    # Manually set up the arrays that collect_cell_matrix would return automatically
    w = []
    r = []
    c = []
 
    if (_is_multifield_fe_basis_component(dv))
       nfields=_nfields(dv)
       fieldid=_fieldid(dv)
       U_Ωbg_cut_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),U_Ωbg_cut_cell_dof_ids)
   end
 
    ## (I) Let's set up the div_dv_div_du_mat_contribs 
    # (first part of the integral, that is ∫( (div_dv)*(div_du)*dΩbg_cut_cells
    div_dv_div_du_mat_contribs=get_array(∫(γ*(∇⋅dv)⋅(∇⋅du))*dΩbg_cut_cells)
    push!(w, div_dv_div_du_mat_contribs)
    push!(r, U_Ωbg_cut_cell_dof_ids)
    push!(c, U_Ωbg_cut_cell_dof_ids)
 
    ## (II) Set up  ∫( (div_dv)*(-div_du_l2_proj_agg_cells))*dΩbg_cut_cells
    div_du=∇⋅(_get_single_field_fe_basis(du))
 
    # Compute Π_Q_bb(div_du)
    dqbb_Ωagg_cells =change_domain_bb_to_agg_cells(qbb,ref_agg_cell_to_ref_bb_map,Ωbg_agg_cells,agg_cells_to_aggregate)  
    agg_cells_rhs_contribs=get_array(∫(dqbb_Ωagg_cells⋅div_du)dΩbg_agg_cells) 
    ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(U_agg_cells_local_dof_ids,agg_cells_rhs_contribs)
    rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)
    div_dv_l2_proj_bb_dofs=lazy_map(\,p_lhs,rhs)
 
    # Generate bb-wise array of fields. For each aggregate's bounding box,
    # it provides the l2 projection of all basis functions in dp 
    # restricted to the cells included in the bounding box of the aggregate  
    div_dv_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination, div_dv_l2_proj_bb_dofs, Gridap.CellData.get_data(qbb))
 
    div_dv_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),lazy_map(Reindex(div_dv_l2_proj_bb_array),agg_cells_to_aggregate),ref_agg_cell_to_ref_bb_map)
    div_du_l2_proj_bb_array_agg_cells=lazy_map(transpose, div_dv_l2_proj_bb_array_agg_cells)
 
    if (_is_multifield_fe_basis_component(du))
        @assert _is_multifield_fe_basis_component(du)
        @assert _nfields(du)==_nfields(dv)
        nfields=_nfields(du)
        fieldid=_fieldid(dv)
        div_du_l2_proj_bb_array_agg_cells=lazy_map(Gridap.Fields.BlockMap((fieldid,nfields),fieldid),div_du_l2_proj_bb_array_agg_cells)
    end 
 
    div_du_l2_proj_agg_cells = Gridap.CellData.GenericCellField(div_du_l2_proj_bb_array_agg_cells,Ωbg_agg_cells,ReferenceDomain())
 
    # In the MultiField case, I have had to add this change domain 
    # call before setting up the term right below. Otherwise, we get an error 
    # when trying to multiply the fields. Not sure why this is happening
    # dv_Ωagg_cells=Gridap.CellData.change_domain(dv,Ωagg_cells,ReferenceDomain())
    div_dv_Ωagg_cells=Gridap.CellData.change_domain(∇⋅dv,Ωbg_agg_cells,ReferenceDomain()) 
 
    proj_div_dv_div_du_mat_contribs=get_array(∫(γ*(-1.0)*div_dv_Ωagg_cells*(div_du_l2_proj_agg_cells))*dΩbg_cut_cells)
 
    if (_is_multifield_fe_basis_component(du))
       nfields=_nfields(du)
       fieldid=_fieldid(du)
       U_cut_cells_to_aggregate_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),U_cut_cells_to_aggregate_dof_ids)
    end
 
    push!(w, proj_div_dv_div_du_mat_contribs)
    push!(r, U_Ωbg_cut_cell_dof_ids)           # test
    push!(c, U_cut_cells_to_aggregate_dof_ids) # trial 
    w,r,c
end

function set_up_bulk_ghost_penalty_lhs(aggregate_to_local_cells,
    agg_cells_to_aggregate,
    ref_agg_cell_to_ref_bb_map,
    dΩagg_cells,
    ubb, 
    vbb)
Ωagg_cells=dΩagg_cells.quad.trian

# Change domain of vbb (test) from Ωbb to Ωagg_cells
vbb_Ωagg_cells=change_domain_bb_to_agg_cells(vbb,
              ref_agg_cell_to_ref_bb_map,
              Ωagg_cells,
              agg_cells_to_aggregate)

# Change domain of ubb (trial) from Ωbb to Ωagg_cells                                             
ubb_Ωagg_cells=change_domain_bb_to_agg_cells(ubb,
             ref_agg_cell_to_ref_bb_map,
             Ωagg_cells,
             agg_cells_to_aggregate)

# Compute contributions to LHS of L2 projection
agg_cells_to_lhs_contribs=get_array(∫(vbb_Ωagg_cells⋅ubb_Ωagg_cells)dΩagg_cells)

# Finally assemble LHS contributions
ass_lhs_map=BulkGhostPenaltyAssembleLhsMap(agg_cells_to_lhs_contribs)
lazy_map(ass_lhs_map,aggregate_to_local_cells)
end 

"""
  dv, du: Test and trial basis functions. They may be components of a MultiFieldCellField

  # Compute and assemble the bulk penalty stabilization term
  # ∫( (div_dv-div_dv_l2_proj_agg_cells)*(div_du-div_du_l2_proj_agg_cells))*dΩ_cut_cells
"""
function div_penalty_stabilization_collect_cell_matrix_on_cut_cells_full(aggregate_to_local_cells,
                                                              agg_cells_to_aggregate,
                                                              ref_agg_cell_to_ref_bb_map,
                                                              dΩbg_agg_cells,
                                                              dv,    # Test basis
                                                              du,    # Trial basis (to project)
                                                              qbb,  # Bounding box space test basis
                                                              p_lhs,
                                                              U_Ωbg_cut_cell_dof_ids,
                                                              U_agg_cells_local_dof_ids,
                                                              U_cut_cells_to_aggregate_dof_ids,
                                                              γ,
                                                              dΩbg_cut_cells) 

    Ωbg_agg_cells=dΩbg_agg_cells.quad.trian

    ## (I) Compute projections div_dv_l2_proj_agg_cells & div_du_l2_proj_agg_cells

    div_du=∇⋅(_get_single_field_fe_basis(du))

    # Change domain of vbb (test) from Ωbb to Ωagg_cells
    dqbb_Ωagg_cells =change_domain_bb_to_agg_cells(qbb,ref_agg_cell_to_ref_bb_map,Ωbg_agg_cells,agg_cells_to_aggregate)  
    agg_cells_rhs_contribs=get_array(∫(dqbb_Ωagg_cells⋅div_du)dΩbg_agg_cells) 
    ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(U_agg_cells_local_dof_ids,agg_cells_rhs_contribs)
    rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)

    # TO-DO: optimize using our own optimized version Gridap.Fields.Map 
    # of backslash that re-uses storage for lu factors among cells, etc.
    div_dv_l2_proj_bb_dofs=lazy_map(\,p_lhs,rhs)


    # Generate bb-wise array of fields. For each aggregate's bounding box,
    # it provides the l2 projection of all basis functions in du 
    # restricted to the cells included in the bounding box of the aggregate   
    div_dv_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination, div_dv_l2_proj_bb_dofs, Gridap.CellData.get_data(qbb))

    # # Change domain of dv_l2_proj_bb_array from bb to agg_cells
    div_dv_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),lazy_map(Reindex(div_dv_l2_proj_bb_array),agg_cells_to_aggregate),ref_agg_cell_to_ref_bb_map)
    div_du_l2_proj_bb_array_agg_cells=lazy_map(transpose, div_dv_l2_proj_bb_array_agg_cells)
    
    if (_is_multifield_fe_basis_component(du))
      @assert _is_multifield_fe_basis_component(du)
      @assert _nfields(du)==_nfields(dv)
      nfields=_nfields(du)
      fieldid=_fieldid(dv)
      div_du_l2_proj_bb_array_agg_cells=lazy_map(Gridap.Fields.BlockMap((fieldid,nfields),fieldid),div_du_l2_proj_bb_array_agg_cells)
    end 

    if (_is_multifield_fe_basis_component(dv))
        @assert _is_multifield_fe_basis_component(dv)
        @assert _nfields(du)==_nfields(dv)
        nfields=_nfields(dv)
        fieldid=_fieldid(dv)
        div_dv_l2_proj_bb_array_agg_cells=lazy_map(
                                    Gridap.Fields.BlockMap(nfields,fieldid),
                                    div_dv_l2_proj_bb_array_agg_cells)
    end 

    div_du_l2_proj_agg_cells = Gridap.CellData.GenericCellField(div_du_l2_proj_bb_array_agg_cells,Ωbg_agg_cells,ReferenceDomain())
    div_dv_l2_proj_agg_cells = Gridap.CellData.GenericCellField(div_dv_l2_proj_bb_array_agg_cells,Ωbg_agg_cells,ReferenceDomain())

   #  # BlockMap preparatory steps for the dof ids
    if (_is_multifield_fe_basis_component(dv))
         nfields=_nfields(dv)
         fieldid=_fieldid(dv)
         U_Ωbg_cut_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),U_Ωbg_cut_cell_dof_ids)
    end

    if (_is_multifield_fe_basis_component(du))
         nfields=_nfields(du)
         fieldid=_fieldid(du)
         U_cut_cells_to_aggregate_dof_ids=
         lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),U_cut_cells_to_aggregate_dof_ids)
    end

    # Manually set up the arrays that collect_cell_matrix would return automatically
    w = []
    r = []
    c = []

    ## (1) div_dv*div_du term
    div_dv_div_du_mat_contribs=get_array(∫(γ*(∇⋅dv)⋅(∇⋅du))*dΩbg_cut_cells)
    push!(w, div_dv_div_du_mat_contribs)
    push!(r, U_Ωbg_cut_cell_dof_ids) 
    push!(c, U_Ωbg_cut_cell_dof_ids)

    ## (2) div_dv*proj_div_du term
    div_dv_Ωagg_cells=Gridap.CellData.change_domain(∇⋅dv,Ωbg_agg_cells,ReferenceDomain()) 
    proj_div_dv_div_du_mat_contribs=get_array(∫(γ*(-1.0)*div_dv_Ωagg_cells⋅(div_du_l2_proj_agg_cells))*dΩbg_cut_cells)

    push!(w, proj_div_dv_div_du_mat_contribs)
    push!(r, U_Ωbg_cut_cell_dof_ids) 
    push!(c, U_cut_cells_to_aggregate_dof_ids)

    ## (3) proj_div_dv*proj_div_du
    proj_div_dv_proj_div_du_mat_contribs=get_array(∫(γ*div_dv_l2_proj_agg_cells⋅(div_du_l2_proj_agg_cells))*dΩbg_cut_cells)

    push!(w, proj_div_dv_proj_div_du_mat_contribs)
    push!(r, U_cut_cells_to_aggregate_dof_ids) 
    push!(c, U_cut_cells_to_aggregate_dof_ids)

    ## (4) proj_div_dv*div_du
    div_du_Ωagg_cells=Gridap.CellData.change_domain(∇⋅du,Ωbg_agg_cells,ReferenceDomain()) 
    proj_div_dv_div_du_mat_contribs=get_array(∫(-γ*div_dv_l2_proj_agg_cells⋅div_du_Ωagg_cells)*dΩbg_cut_cells)

    push!(w, proj_div_dv_div_du_mat_contribs)
    push!(r, U_cut_cells_to_aggregate_dof_ids) 
    push!(c, U_Ωbg_cut_cell_dof_ids)

    w, r, c
end

## Create collect function on cut cells only
"""
  dv, du: Test and trial basis functions. They may be components of a MultiFieldCellField

  # Compute and assemble the bulk penalty stabilization term
  # ∫( (div_dv)*(dp-dp_l2_proj_agg_cells))*dΩ_cut_cells (reduced)
    ##### ∫( (div_dv-div_dv_l2_proj_agg_cells)*(dp-dp_l2_proj_agg_cells))*dΩ_cut_cells (full form, not implented here)

"""
function pmix_penalty_stabilization_collect_cell_matrix_on_cut_cells(aggregate_to_local_cells,
                                                              agg_cells_to_aggregate,
                                                              ref_agg_cell_to_ref_bb_map,
                                                              dΩbg_agg_cells,
                                                              dq,    # Test basis
                                                              dp,    # Trial basis (to project)
                                                              dqbb,  # Bounding box space test basis
                                                              p_lhs,
                                                              U_Ωbg_cut_cell_dof_ids,
                                                              P_Ωbg_cut_cell_dof_ids,
                                                              P_agg_cells_local_dof_ids,
                                                              P_cut_cells_to_aggregate_dof_ids,
                                                              γ,
                                                              dΩbg_cut_cells,
                                                              dv) 


    Ωbg_agg_cells=dΩbg_agg_cells.quad.trian

    ## (I) Compute projections dq_l2_proj_agg_cells & dp_l2_proj_agg_cells
    # Change domain of qbb (test) from Ωbb to Ωagg_cells
    dqbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(dqbb,
                                                  ref_agg_cell_to_ref_bb_map,
                                                  Ωbg_agg_cells,
                                                  agg_cells_to_aggregate)             

    dp_single_field=_get_single_field_fe_basis(dp)                                              
    agg_cells_rhs_contribs=get_array(∫(dqbb_Ωbg_agg_cells⋅dp_single_field)dΩbg_agg_cells)
    ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(P_agg_cells_local_dof_ids,agg_cells_rhs_contribs)
    rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)

    # TO-DO: optimize using our own optimized version Gridap.Fields.Map 
    # of backslash that re-uses storage for lu factors among cells, etc.
    dq_l2_proj_bb_dofs=lazy_map(\,p_lhs,rhs)

    # Generate bb-wise array of fields. For each aggregate's bounding box,
    # it provides the l2 projection of all basis functions in dp 
    # restricted to the cells included in the bounding box of the aggregate   
    dq_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination,
                                dq_l2_proj_bb_dofs,
                                Gridap.CellData.get_data(dqbb))

    # # Change domain of dq_l2_proj_bb_array from bb to agg_cells
    dq_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),
                                           lazy_map(Reindex(dq_l2_proj_bb_array),agg_cells_to_aggregate),
                                           ref_agg_cell_to_ref_bb_map)

    dp_l2_proj_bb_array_agg_cells=lazy_map(transpose, dq_l2_proj_bb_array_agg_cells)
    
    if (_is_multifield_fe_basis_component(dp))
        @assert _is_multifield_fe_basis_component(dq)
        @assert _nfields(dp)==_nfields(dq)
        nfields=_nfields(dp)
        fieldid=_fieldid(dp)
        dp_l2_proj_bb_array_agg_cells=lazy_map(
                                    Gridap.Fields.BlockMap((1,nfields),fieldid),
                                    dp_l2_proj_bb_array_agg_cells)
    end

    # if (_is_multifield_fe_basis_component(dq))
    #     @assert _is_multifield_fe_basis_component(dq)
    #     @assert _nfields(dp)==_nfields(dq)
    #     nfields=_nfields(dq)
    #     fieldid=_fieldid(dq)
    #     dq_l2_proj_bb_array_agg_cells=lazy_map(
    #                                 Gridap.Fields.BlockMap(nfields,fieldid),
    #                                 dq_l2_proj_bb_array_agg_cells)
    #  end 

    dp_l2_proj_agg_cells = Gridap.CellData.GenericCellField(dp_l2_proj_bb_array_agg_cells, Ωbg_agg_cells,ReferenceDomain()) 
    # dq_l2_proj_agg_cells = Gridap.CellData.GenericCellField(dq_l2_proj_bb_array_agg_cells, Ωbg_agg_cells,ReferenceDomain()) 

    # BlockMap preparatory steps for the dof ids
    if (_is_multifield_fe_basis_component(dq))
        nfields=_nfields(dq)
        fieldid=_fieldid(dq)
        P_Ωbg_cut_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),P_Ωbg_cut_cell_dof_ids)
    end

    if (_is_multifield_fe_basis_component(dp))
        nfields=_nfields(dp)
        fieldid=_fieldid(dp)
        P_cut_cells_to_aggregate_dof_ids=
        lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),P_cut_cells_to_aggregate_dof_ids)
    end

    if (_is_multifield_fe_basis_component(dv))
        nfields=_nfields(dv)
        fieldid=_fieldid(dv)
        U_Ωbg_cut_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),U_Ωbg_cut_cell_dof_ids)
    end

    # Manually set up the arrays that collect_cell_matrix would return automatically
    w = []
    r = []
    c = []

    ## (1) div_dv*dp term
    # div_dv_Ωagg_cells=Gridap.CellData.change_domain(∇⋅dv,Ωbg_agg_cells,ReferenceDomain()) 
    # div_dv_dp_mat_contribs=get_array(∫(γ*div_dv_Ωagg_cells⋅dp)*dΩbg_cut_cells)
    div_dv_dp_mat_contribs=get_array(∫(γ*(∇⋅dv)⋅dp)*dΩbg_cut_cells)

    push!(w, div_dv_dp_mat_contribs)
    push!(r, U_Ωbg_cut_cell_dof_ids) 
    push!(c, P_Ωbg_cut_cell_dof_ids)

    ## (2) div_dv*proj_dp term
    div_dv_Ωagg_cells=Gridap.CellData.change_domain(∇⋅dv,Ωbg_agg_cells,ReferenceDomain()) 
    div_dv_proj_dp_mat_contribs=get_array(∫(γ*(-1.0)*div_dv_Ωagg_cells⋅(dp_l2_proj_agg_cells))*dΩbg_cut_cells)

    push!(w, div_dv_proj_dp_mat_contribs)    
    push!(r, U_Ωbg_cut_cell_dof_ids) 
    push!(c, P_cut_cells_to_aggregate_dof_ids)

    # ## (3) proj_dq*proj_dp
    # proj_dq_proj_dp_mat_contribs=get_array(∫(γ*dq_l2_proj_agg_cells⋅(dp_l2_proj_agg_cells))*dΩbg_cut_cells)

    # push!(w, proj_dq_proj_dp_mat_contribs)
    # push!(r, cut_cells_to_aggregate_dof_ids) 
    # push!(c, cut_cells_to_aggregate_dof_ids)

    # ## (4) proj_dq*dp
    # dpΩagg_cells=Gridap.CellData.change_domain(dp,Ωbg_agg_cells,ReferenceDomain())
    # proj_dq_dp_mat_contribs=get_array(∫(-γ*dq_l2_proj_agg_cells⋅dpΩagg_cells)*dΩbg_cut_cells)

    # push!(w, proj_dq_dp_mat_contribs)
    # push!(r, cut_cells_to_aggregate_dof_ids) 
    # push!(c, Ωbg_cut_cell_dof_ids)

    w, r, c
end

## Create collect function on cut cells only
"""
  dv, du: Test and trial basis functions. They may be components of a MultiFieldCellField

  # Compute and assemble the bulk penalty stabilization term
  # ∫( (div_dv)*(dp-dp_l2_proj_agg_cells))*dΩ_cut_cells (reduced)
    ##### ∫( (div_dv-div_dv_l2_proj_agg_cells)*(dp-dp_l2_proj_agg_cells))*dΩ_cut_cells (full form, not implented here)

"""
function dmix_penalty_stabilization_collect_cell_matrix_on_cut_cells(aggregate_to_local_cells,
                                                              agg_cells_to_aggregate,
                                                              ref_agg_cell_to_ref_bb_map,
                                                              dΩbg_agg_cells,
                                                              dq,    # Test basis
                                                              dp,    # Trial basis (to project)
                                                              dqbb,  # Bounding box space test basis
                                                              p_lhs,
                                                              U_Ωbg_cut_cell_dof_ids,
                                                              P_Ωbg_cut_cell_dof_ids,
                                                              P_agg_cells_local_dof_ids,
                                                              P_cut_cells_to_aggregate_dof_ids,
                                                              γ,
                                                              dΩbg_cut_cells,
                                                              du) 


    Ωbg_agg_cells=dΩbg_agg_cells.quad.trian

    ## (I) Compute projections dq_l2_proj_agg_cells & dp_l2_proj_agg_cells
    # Change domain of qbb (test) from Ωbb to Ωagg_cells
    dqbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(dqbb,
                                                  ref_agg_cell_to_ref_bb_map,
                                                  Ωbg_agg_cells,
                                                  agg_cells_to_aggregate)             

    dp_single_field=_get_single_field_fe_basis(dp)                                              
    agg_cells_rhs_contribs=get_array(∫(dqbb_Ωbg_agg_cells⋅dp_single_field)dΩbg_agg_cells)
    ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(P_agg_cells_local_dof_ids,agg_cells_rhs_contribs)
    rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)

    # TO-DO: optimize using our own optimized version Gridap.Fields.Map 
    # of backslash that re-uses storage for lu factors among cells, etc.
    dq_l2_proj_bb_dofs=lazy_map(\,p_lhs,rhs)

    # Generate bb-wise array of fields. For each aggregate's bounding box,
    # it provides the l2 projection of all basis functions in dp 
    # restricted to the cells included in the bounding box of the aggregate   
    dq_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination,
                                dq_l2_proj_bb_dofs,
                                Gridap.CellData.get_data(dqbb))

    # # Change domain of dq_l2_proj_bb_array from bb to agg_cells
    dq_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),
                                           lazy_map(Reindex(dq_l2_proj_bb_array),agg_cells_to_aggregate),
                                           ref_agg_cell_to_ref_bb_map)

    if (_is_multifield_fe_basis_component(dq))
        @assert _is_multifield_fe_basis_component(dq)
        @assert _nfields(dp)==_nfields(dq)
        nfields=_nfields(dq)
        fieldid=_fieldid(dq)
        dq_l2_proj_bb_array_agg_cells=lazy_map(
                                    Gridap.Fields.BlockMap(nfields,fieldid),
                                    dq_l2_proj_bb_array_agg_cells)
     end 

    # dp_l2_proj_agg_cells = Gridap.CellData.GenericCellField(dp_l2_proj_bb_array_agg_cells, Ωbg_agg_cells,ReferenceDomain()) 
    dq_l2_proj_agg_cells = Gridap.CellData.GenericCellField(dq_l2_proj_bb_array_agg_cells, Ωbg_agg_cells,ReferenceDomain()) 

    # BlockMap preparatory steps for the dof ids
    if (_is_multifield_fe_basis_component(dq))
        nfields=_nfields(dq)
        fieldid=_fieldid(dq)
        P_Ωbg_cut_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),P_Ωbg_cut_cell_dof_ids)
    end

    if (_is_multifield_fe_basis_component(dp))
        nfields=_nfields(dp)
        fieldid=_fieldid(dp)
        P_cut_cells_to_aggregate_dof_ids=
        lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),P_cut_cells_to_aggregate_dof_ids)
    end

    if (_is_multifield_fe_basis_component(du))
        nfields=_nfields(du)
        fieldid=_fieldid(du)
        U_Ωbg_cut_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),U_Ωbg_cut_cell_dof_ids)
    end

    # Manually set up the arrays that collect_cell_matrix would return automatically
    w = []
    r = []
    c = []

    ## (1) div_du*dq term
    # div_dv_Ωagg_cells=Gridap.CellData.change_domain(∇⋅dv,Ωbg_agg_cells,ReferenceDomain()) 
    # div_dv_dp_mat_contribs=get_array(∫(γ*div_dv_Ωagg_cells⋅dp)*dΩbg_cut_cells)
    div_dv_dp_mat_contribs=get_array(∫(γ*(∇⋅du)⋅dq)*dΩbg_cut_cells)

    push!(w, div_dv_dp_mat_contribs)
    push!(r, P_Ωbg_cut_cell_dof_ids) 
    push!(c, U_Ωbg_cut_cell_dof_ids)

    ## (2) div_du*proj_dq term
    div_du_Ωagg_cells=Gridap.CellData.change_domain(∇⋅du,Ωbg_agg_cells,ReferenceDomain()) 
    div_dv_proj_dp_mat_contribs=get_array(∫(γ*(-1.0)*div_du_Ωagg_cells⋅(dq_l2_proj_agg_cells))*dΩbg_cut_cells)

    push!(w, div_dv_proj_dp_mat_contribs)    
    push!(r, P_cut_cells_to_aggregate_dof_ids)
    push!(c, U_Ωbg_cut_cell_dof_ids) 


    # ## (3) proj_dq*proj_dp
    # proj_dq_proj_dp_mat_contribs=get_array(∫(γ*dq_l2_proj_agg_cells⋅(dp_l2_proj_agg_cells))*dΩbg_cut_cells)

    # push!(w, proj_dq_proj_dp_mat_contribs)
    # push!(r, cut_cells_to_aggregate_dof_ids) 
    # push!(c, cut_cells_to_aggregate_dof_ids)

    # ## (4) proj_dq*dp
    # dpΩagg_cells=Gridap.CellData.change_domain(dp,Ωbg_agg_cells,ReferenceDomain())
    # proj_dq_dp_mat_contribs=get_array(∫(-γ*dq_l2_proj_agg_cells⋅dpΩagg_cells)*dΩbg_cut_cells)

    # push!(w, proj_dq_dp_mat_contribs)
    # push!(r, cut_cells_to_aggregate_dof_ids) 
    # push!(c, Ωbg_cut_cell_dof_ids)

    w, r, c
end

## Create collect function on cut cells only for rhs 
## TODO: rhs_g_func does not seem to work when function is passed on
"""
  dv, du: Test and trial basis functions. They may be components of a MultiFieldCellField

  # Compute and assemble the bulk penalty stabilization term
  # ∫( (div_dv)*(dp-dp_l2_proj_agg_cells))*dΩ_cut_cells (reduced)
    ##### ∫( (div_dv-div_dv_l2_proj_agg_cells)*(dp-dp_l2_proj_agg_cells))*dΩ_cut_cells (full form, not implented here)

"""
function dmix_penalty_stabilization_collect_cell_vector_on_cut_cells(aggregate_to_local_cells,
                                                              agg_cells_to_aggregate,
                                                              ref_agg_cell_to_ref_bb_map,
                                                              dΩbg_agg_cells,
                                                              dq,    # Test basis
                                                              dp,    # Trial basis (to project)
                                                              dqbb,  # Bounding box space test basis
                                                              p_lhs,
                                                              P_Ωbg_cut_cell_dof_ids,
                                                              P_agg_cells_local_dof_ids,
                                                              P_cut_cells_to_aggregate_dof_ids,
                                                              γ,
                                                              dΩbg_cut_cells,
                                                              rhs_g_func) 


    Ωbg_agg_cells=dΩbg_agg_cells.quad.trian

    ## Compute projections dq_l2_proj_agg_cells & dp_l2_proj_agg_cells
    # Change domain of qbb (test) from Ωbb to Ωagg_cells
    dqbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(dqbb,
                                                  ref_agg_cell_to_ref_bb_map,
                                                  Ωbg_agg_cells,
                                                  agg_cells_to_aggregate)             

    dp_single_field=_get_single_field_fe_basis(dp)                                              
    agg_cells_rhs_contribs=get_array(∫(dqbb_Ωbg_agg_cells⋅dp_single_field)dΩbg_agg_cells)
    ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(P_agg_cells_local_dof_ids,agg_cells_rhs_contribs)
    rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)

    # TO-DO: optimize using our own optimized version Gridap.Fields.Map 
    # of backslash that re-uses storage for lu factors among cells, etc.
    dq_l2_proj_bb_dofs=lazy_map(\,p_lhs,rhs)

    # Generate bb-wise array of fields. For each aggregate's bounding box,
    # it provides the l2 projection of all basis functions in dp 
    # restricted to the cells included in the bounding box of the aggregate   
    dq_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination,
                                dq_l2_proj_bb_dofs,
                                Gridap.CellData.get_data(dqbb))

    # # Change domain of dq_l2_proj_bb_array from bb to agg_cells
    dq_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),
                                           lazy_map(Reindex(dq_l2_proj_bb_array),agg_cells_to_aggregate),
                                           ref_agg_cell_to_ref_bb_map)

    if (_is_multifield_fe_basis_component(dq))
        @assert _is_multifield_fe_basis_component(dq)
        @assert _nfields(dp)==_nfields(dq)
        nfields=_nfields(dq)
        fieldid=_fieldid(dq)
        dq_l2_proj_bb_array_agg_cells=lazy_map(
                                    Gridap.Fields.BlockMap(nfields,fieldid),
                                    dq_l2_proj_bb_array_agg_cells)
     end 

    dq_l2_proj_agg_cells = Gridap.CellData.GenericCellField(dq_l2_proj_bb_array_agg_cells, Ωbg_agg_cells,ReferenceDomain()) 

    # BlockMap preparatory steps for the dof ids
    if (_is_multifield_fe_basis_component(dq))
        nfields=_nfields(dq)
        fieldid=_fieldid(dq)
        P_Ωbg_cut_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),P_Ωbg_cut_cell_dof_ids)
    end

    if (_is_multifield_fe_basis_component(dq))
        nfields=_nfields(dq)
        fieldid=_fieldid(dq)
        P_cut_cells_to_aggregate_dof_ids=
        lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),P_cut_cells_to_aggregate_dof_ids)
    end

    # Manually set up the arrays that collect_cell_vector would return automatically
    w = []
    r = []

    ## (1) rhs_g*dq term
    rhs_g_dp_mat_contribs=get_array(∫((rhs_g_func⋅dq)*γ)*dΩbg_cut_cells)

    push!(w, rhs_g_dp_mat_contribs)
    push!(r, P_Ωbg_cut_cell_dof_ids) 

    # (2) div_du*proj_dq term
    rhs_g_proj_dp_mat_contribs=get_array(∫((rhs_g_func⋅dq_l2_proj_agg_cells)*γ*(-1.0))*dΩbg_cut_cells)

    push!(w, rhs_g_proj_dp_mat_contribs)    
    push!(r, P_cut_cells_to_aggregate_dof_ids)

    w, r
end