function assemble_darcy_preconditioner_matrix(cutgeo, degree, X, Y, s, h, β₀,
                                              aggregate_to_local_cells,
                                              agg_cells_to_aggregate, 
                                              ref_agg_cell_to_ref_bb_map,
                                              dΩbg_agg_cells,
                                              dΩbg_cut_cells,
                                              # FLUX-RELATED DATA
                                              dv,    # Test basis
                                              du,    # Trial basis (to project)
                                              vbb,  # Bounding box space test basis
                                              u_lhs,
                                              U_Ωbg_cut_cell_dof_ids, 
                                              U_agg_cells_local_dof_ids,
                                              U_cut_cells_to_aggregate_dof_ids,
                                              γ_u,
                                              # PRESSURE-RELATED DATA
                                              dq,    # Test basis
                                              dp,    # Trial basis (to project)
                                              qbb,  # Bounding box space test basis
                                              p_lhs,
                                              P_Ωbg_cut_cell_dof_ids, 
                                              P_agg_cells_local_dof_ids,
                                              P_cut_cells_to_aggregate_dof_ids,
                                              γ_p,

                                              γ_div)

   
   # Physical domain
   Ω = Triangulation(cutgeo,PHYSICAL)
   dΩ = Measure(Ω,degree)
   
   # Bounadary
   ΓN    = EmbeddedBoundary(cutgeo)
   dΓN   = Measure(ΓN,degree)
   nN    = get_normal_vector(ΓN)
   β     = β₀*h.^(s)
   a((u,p), (v,q)) = ∫(v⋅u)dΩ + ∫(p⋅q)dΩ + ∫(β*(u⋅nN)*(v⋅nN))dΓN    
   
   dx = get_trial_fe_basis(X)
   dy = get_fe_basis(Y)

   
   ## FULL stabilization terms:
   wu_full,ru_full,cu_full=
     bulk_ghost_penalty_stabilization_collect_cell_matrix_on_cut_cells_full(aggregate_to_local_cells,
                                                        agg_cells_to_aggregate,
                                                        ref_agg_cell_to_ref_bb_map,
                                                        dΩbg_agg_cells,
                                                        dv,    # Test basis
                                                        du,    # Trial basis (to project)
                                                        vbb,  # Bounding box space test basis
                                                        u_lhs,
                                                        U_Ωbg_cut_cell_dof_ids, 
                                                        U_agg_cells_local_dof_ids,
                                                        U_cut_cells_to_aggregate_dof_ids,
                                                        γ_u,
                                                        dΩbg_cut_cells)                           

   wp_full,rp_full,cp_full=
      bulk_ghost_penalty_stabilization_collect_cell_matrix_on_cut_cells_full(aggregate_to_local_cells,
                                                        agg_cells_to_aggregate,
                                                        ref_agg_cell_to_ref_bb_map,
                                                        dΩbg_agg_cells,
                                                        dq,    # Test basis
                                                        dp,    # Trial basis (to project)
                                                        qbb,  # Bounding box space test basis
                                                        p_lhs,
                                                        P_Ωbg_cut_cell_dof_ids, 
                                                        P_agg_cells_local_dof_ids,
                                                        P_cut_cells_to_aggregate_dof_ids,
                                                        γ_p,
                                                        dΩbg_cut_cells)   

   wdiv_full, rdiv_full, cdiv_full = 
      div_penalty_stabilization_collect_cell_matrix_on_cut_cells_full(aggregate_to_local_cells,
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
                                                        γ_div,
                                                        dΩbg_cut_cells)

     assem=SparseMatrixAssembler(X,Y)
     wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
     push!(wrc[1], wu_full...)
     push!(wrc[2], ru_full...)
     push!(wrc[3], cu_full...)

     push!(wrc[1], wp_full...)
     push!(wrc[2], rp_full...)
     push!(wrc[3], cp_full...)

     push!(wrc[1], wdiv_full...)
     push!(wrc[2], rdiv_full...)
     push!(wrc[3], cdiv_full...)

     assemble_matrix(assem, wrc)
end