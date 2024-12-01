using Gridap
using GridapEmbedded
using FillArrays
using LinearAlgebra

include("aggregates_bounding_boxes_tools.jl")
include("bulk_ghost_penalty_stab_tools.jl")
include("fields_and_blocks_tools.jl")
include("BulkGhostPenaltyAssembleMaps.jl")

# Manufactured solution
order = 0
uex(x) = -VectorValue(2*x[1],2*x[2])
pex(x) = (x[1]^2 + x[2]^2)
divuex(x) = -4.0

# Select geometry
R = 0.2
geom = disk(R, x0=Point(0.5,0.5))

# Setup background model
n=16
partition = (n,n)
box = get_metadata(geom)
bgmodel = CartesianDiscreteModel((0,1,0,1),partition)
dp = box.pmax - box.pmin
h = dp[1]/n

# Cut the background model with the mesh
cutgeo = cut(bgmodel,geom)

# Compute mapping among background model 
# cut cells and interior cells 
strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutgeo)

aggregate_to_cells=setup_aggregate_to_cells(aggregates)
aggregates_bounding_box_model=
       setup_aggregates_bounding_box_model(bgmodel,aggregate_to_cells)

# Triangulations
Ωbg  = Triangulation(bgmodel)
Ωact = Triangulation(cutgeo,ACTIVE)

# Physical domain
Ω = Triangulation(cutgeo,PHYSICAL)
degree=2*(order+1)
dΩ = Measure(Ω,degree)

# Set up global spaces 
V = FESpace(Ωact, ReferenceFE(raviart_thomas,Float64,order),conformity=:HDiv)
Q = FESpace(Ωact, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
U = TrialFESpace(V)
P = TrialFESpace(Q)
Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])
dx    = get_trial_fe_basis(X)
dy    = get_fe_basis(Y)
du,dp = dx  
dv,dq = dy 

# Aggregate cells and triangulation
agg_cells    =setup_agg_cells(aggregate_to_cells)
Ωbg_agg_cells=view(Ωbg,agg_cells)

# ref_agg_cell_to_agg_cell_map: \hat{K} -> K 
ref_agg_cell_to_agg_cell_map=get_cell_map(Ωbg_agg_cells)
agg_cells_to_aggregate      =setup_agg_cells_to_aggregate(aggregate_to_cells)
ref_agg_cell_to_ref_bb_map  =setup_ref_agg_cell_to_ref_bb_map(aggregates_bounding_box_model,
                                                            agg_cells_to_aggregate)

# Spaces on bounding boxes                                                            
reffeₚ_bb =ReferenceFE(lagrangian,Float64,order) # Here we MUST use a Q space (not a P space!)
Qbb=FESpace(aggregates_bounding_box_model,reffeₚ_bb,conformity=:L2) # We need a DG space to represent the L2 projection
Pbb=TrialFESpace(Qbb)
pbb=get_trial_fe_basis(Pbb)
qbb=get_fe_basis(Qbb)
reffeᵤ_bb=ReferenceFE(raviart_thomas,Float64,order)
Vbb=FESpace(aggregates_bounding_box_model,reffeᵤ_bb,conformity=:L2)
Ubb=TrialFESpace(Vbb)
ubb=get_trial_fe_basis(Ubb)
vbb=get_fe_basis(Vbb)

# Numerical integration (Measures)
degree=2*(order+1)
dΩbg_agg_cells = Measure(Ωbg_agg_cells,degree)

# LHS of L2 projection on bounding boxes.
aggregate_to_local_cells=setup_aggregate_to_local_cells(aggregate_to_cells)
p_lhs=set_up_bulk_ghost_penalty_lhs(aggregate_to_local_cells,
                                  agg_cells_to_aggregate,
                                  ref_agg_cell_to_ref_bb_map,
                                  dΩbg_agg_cells,
                                  pbb,
                                  qbb)

u_lhs=set_up_bulk_ghost_penalty_lhs(aggregate_to_local_cells,
                                  agg_cells_to_aggregate,
                                  ref_agg_cell_to_ref_bb_map,
                                  dΩbg_agg_cells,
                                  ubb,
                                  vbb)   

# Selecting relevant global dofs ids of aggregate cells (from background mesh)
Ωbg_agg_cell_dof_ids   = get_cell_dof_ids(X,Ωbg_agg_cells)
U_Ωbg_agg_cell_dof_ids = _restrict_to_block(Ωbg_agg_cell_dof_ids, 1) 
P_Ωbg_agg_cell_dof_ids = _restrict_to_block(Ωbg_agg_cell_dof_ids, 2)

# Computing local (per aggregate) dof ids 
U_agg_cells_local_dof_ids=
   compute_agg_cells_local_dof_ids(U_Ωbg_agg_cell_dof_ids, aggregate_to_local_cells)
P_agg_cells_local_dof_ids=
   compute_agg_cells_local_dof_ids(P_Ωbg_agg_cell_dof_ids, aggregate_to_local_cells)

# Compute global dofs ids per aggregate and reindex these 
U_aggregate_dof_ids=compute_aggregate_dof_ids(U_Ωbg_agg_cell_dof_ids,aggregate_to_cells)
U_agg_cells_to_aggregate_dof_ids=lazy_map(Reindex(U_aggregate_dof_ids),agg_cells_to_aggregate)
P_aggregate_dof_ids=compute_aggregate_dof_ids(P_Ωbg_agg_cell_dof_ids,aggregate_to_cells)
P_agg_cells_to_aggregate_dof_ids=lazy_map(Reindex(P_aggregate_dof_ids),agg_cells_to_aggregate)
P_Ωbg_agg_cell_dof_ids

# parameters
γ = 10.0 # Interior bulk-penalty stabilization parameter
h_U = 1.0

###########################################
### STABILIZATION ON Ωagg\Troot ###
###########################################

## (1) FIND ALL INTERIOR CELLS WHICH ARE NOT PART OF AN AGGREGATE
int_nonagg_cell_to_cells = setup_int_nonagg_cell_to_cells(aggregates)
int_nonagg_cells = setup_int_nonagg_cells(int_nonagg_cell_to_cells)

## (2) FIND THE INTERIOR CELLS
Ω_agg_cells = view(Ω,agg_cells)            # cells in aggregates
int_cells   = Ω_agg_cells.parent.b.tface_to_mface 

## (3) FIND THE INTERIOR AGGREGATE (/ROOT) CELLS AND CUT CELLS 
root_cells = setup_root_cells(int_cells, int_nonagg_cells)
cut_cells  = setup_cut_cells(agg_cells, root_cells)

## (4) CUT CELLS AND MEASURE
Ωbg_cut_cells   = view(Ωbg,cut_cells)            # cut cells (part of aggregate)
dΩbg_cut_cells  = Measure(Ωbg_cut_cells,degree)

# (5) DOF IDS related to cut cells
Ωbg_cut_cell_dof_ids   = get_cell_dof_ids(X,Ωbg_cut_cells)
U_Ωbg_cut_cell_dof_ids = _restrict_to_block(Ωbg_cut_cell_dof_ids, 1) 
P_Ωbg_cut_cell_dof_ids = _restrict_to_block(Ωbg_cut_cell_dof_ids, 2)

# (6) Defining cut_cells_to_aggregate_dof_ids
aggregate_to_cut_cells = setup_aggregate_to_cut_cells(aggregates, root_cells)
cut_cells_in_agg_cells_to_aggregate = setup_cut_cells_in_agg_cells_to_aggregate(aggregate_to_cut_cells)
U_cut_cells_to_aggregate_dof_ids=lazy_map(Reindex(U_aggregate_dof_ids),cut_cells_in_agg_cells_to_aggregate)
P_cut_cells_to_aggregate_dof_ids=lazy_map(Reindex(P_aggregate_dof_ids),cut_cells_in_agg_cells_to_aggregate)

# (7) Compute stabilization terms for u and p 
wu,ru,cu=bulk_ghost_penalty_stabilization_collect_cell_matrix_on_cut_cells(agg_cells_to_aggregate,
                                                        ref_agg_cell_to_ref_bb_map,
                                                        dΩbg_agg_cells,
                                                        dv,    # Test basis
                                                        du,    # Trial basis (to project)
                                                        vbb,  # Bounding box space test basis
                                                        u_lhs,
                                                        U_Ωbg_cut_cell_dof_ids, 
                                                        U_agg_cells_local_dof_ids,
                                                        U_cut_cells_to_aggregate_dof_ids,
                                                        γ,
                                                        dΩbg_cut_cells)                           

wp,rp,cp=bulk_ghost_penalty_stabilization_collect_cell_matrix_on_cut_cells(agg_cells_to_aggregate,
                                                        ref_agg_cell_to_ref_bb_map,
                                                        dΩbg_agg_cells,
                                                        dq,    # Test basis
                                                        dp,    # Trial basis (to project)
                                                        qbb,  # Bounding box space test basis
                                                        p_lhs,
                                                        P_Ωbg_cut_cell_dof_ids, 
                                                        P_agg_cells_local_dof_ids,
                                                        P_cut_cells_to_aggregate_dof_ids,
                                                        γ,
                                                        dΩbg_cut_cells)                                                       
                                                        
#DIV stabilization part
wdiv, rdiv, cdiv = div_penalty_stabilization_collect_cell_matrix_on_cut_cells(agg_cells_to_aggregate,
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

## FULL stabilization terms:
wu_full,ru_full,cu_full=bulk_ghost_penalty_stabilization_collect_cell_matrix_on_cut_cells_full(agg_cells_to_aggregate,
                                                        ref_agg_cell_to_ref_bb_map,
                                                        dΩbg_agg_cells,
                                                        dv,    # Test basis
                                                        du,    # Trial basis (to project)
                                                        vbb,  # Bounding box space test basis
                                                        u_lhs,
                                                        U_Ωbg_cut_cell_dof_ids, 
                                                        U_agg_cells_local_dof_ids,
                                                        U_cut_cells_to_aggregate_dof_ids,
                                                        γ,
                                                        dΩbg_cut_cells)                           

wp_full,rp_full,cp_full=bulk_ghost_penalty_stabilization_collect_cell_matrix_on_cut_cells_full(agg_cells_to_aggregate,
                                                        ref_agg_cell_to_ref_bb_map,
                                                        dΩbg_agg_cells,
                                                        dq,    # Test basis
                                                        dp,    # Trial basis (to project)
                                                        qbb,  # Bounding box space test basis
# (k=0,n=16) 1.29e7, 54500.0, 8.79e-26, 1.16e-4, 1.98e-28       # U_FULL + P_FULL + DIVU_FULL
                                                        p_lhs,
                                                        P_Ωbg_cut_cell_dof_ids, 
                                                        P_agg_cells_local_dof_ids,
                                                        P_cut_cells_to_aggregate_dof_ids,
                                                        γ,
                                                        dΩbg_cut_cells)   

## TEST WITH MANUFACTURED SOLUTIONS
a((u,p),(v,q))=∫(v⋅u+q*p+(∇⋅u)*(∇⋅v))dΩ 
# wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
assem=SparseMatrixAssembler(X,Y)
l((v,q))=∫(v⋅uex+q*pex+(∇⋅v)*divuex)dΩ
b = assemble_vector(l, Y)

function compute_quantities(A,b,dΩ)
   cond_A =  cond(Array(A))
   norm_A =  norm(A)
   sol_x = A\b
   xh = FEFunction(X, sol_x)
   uh,ph = xh
   euh = uex-uh
   eph = pex-ph
   edivuh = divuex-(∇⋅uh)
   norm_euh = sum(∫(euh⋅euh)*dΩ)
   norm_eph = sum(∫(eph*eph)*dΩ)
   norm_edivuh = sum(∫(edivuh⋅edivuh)*dΩ)
   return round(cond_A,sigdigits=3), round(norm_A,sigdigits=3), round(norm_euh,sigdigits=3), round(norm_eph,sigdigits=3), round(norm_edivuh,sigdigits=3)
end

############## (k=0,n=16) #####################################################################
# (k=0,n=16) 2.31e9, 6040.0,  1.20e-27, 7.71e-5, 1.77e-30       # -- NO STAB --
# (k=0,n=16) 3.58e8, 6040.0,  2.0 e-28, 7.71e-5, 1.93e-30       # U
# (k=0,n=16) 3.58e8, 6050.0,  1.67e-27, 7.71e-5, 4.26e-30       # U_FULL
# (k=0,n=16) 7.59e6, 55500.0, 1.64e-26, 1.46e-4, 4.20e-29       # U      + P      + DIVU
# (k=0,n=16) 7.59e6, 55500.0, 6.17e-26, 1.45e-4, 9.62e-29       # U_FULL + P      + DIVU
# (k=0,n=16) 7.64e6, 54500.0, 8.79e-26, 1.46e-4, 1.98e-28       # U_FULL + P      + DIVU_FULL (*)
# (k=0,n=16) 1.28e7, 55500.0, 1.64e-26, 1.16e-4, 4.20e-29       # U      + P_FULL + DIVU
# (k=0,n=16) 1.28e7, 55500.0, 6.17e-26, 1.16e-4, 9.62e-29       # U_FULL + P_FULL + DIVU
# (k=0,n=16) 1.29e7, 54500.0, 8.79e-26, 1.16e-4, 1.98e-28       # U_FULL + P_FULL + DIVU_FULL
# (k=0,n=16) 7.64e6, 54500.0, 2.75e-26, 1.46e-4, 1.14e-28       # U      + P      + DIVU_FULL (*)
###############################################################################################

############## (k=2,n=16) #####################################################################
# (k=0,n=16) 2.82e29, 9.21e8, 4.0e-12, 7.79e-28, 2.34e-13        # -- NO STAB --
# (k=0,n=16) 2.82e29, 9.21e8, 2.05e-18, 7.79e-28, 1.66e-19       # U
# (k=0,n=16) 2.82e29, 9.21e8, 1.23e-16, 7.79e-28, 1.39e-18       # U_FULL
# (k=0,n=16) 4.72e16, 9.45e9, 3.08e-14, 4.07e-29, 1.76e-17       # U      + P      + DIVU
# (k=0,n=16) 4.72e16, 9.45e9, 2.52e-11, 4.07e-29, 3.17e-16       # U_FULL + P      + DIVU
# (k=0,n=16) 3.00e18, 9.45e9, 3.08e-14, 1.65e-28, 1.76e-17       # U      + P_FULL + DIVU
# (k=0,n=16) 3.00e18, 9.45e9, 2.52e-11, 1.65e-28, 3.17e-16       # U_FULL + P_FULL + DIVU
###############################################################################################

# CONCLUSION (k=0): lowest κ for stab u, p and divu. For u it does not matter whether it is full or not full. P__full would increase the order of κ by 1. 
# (*)NOTE THAT U(_FUL)L + P + DIVU_FULL does not outperform U(_FULL) + P + DIVU, with κ=7.64e6 and κ=7.59e6, respectively.

res_nostab
res_stab_u
res_stab_ufull
res_stab_updivu
res_stab_ufullpdivu
res_stab_ufullpdivufull
res_stab_upfulldivu 
res_stab_ufullpfulldivu 

# NO STAB
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
A = assemble_matrix(assem, wrc)
res_nostab = compute_quantities(A,b,dΩ)
# (k=0,n=16) 2.31e9, 6040.0, 1.20e-27, 7.71e-5, 1.77e-30        # NO STAB

# ONLY U
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu...)
push!(wrc[2], ru...)
push!(wrc[3], cu...)
A = assemble_matrix(assem, wrc)
res_stab_u = compute_quantities(A,b,dΩ)
# (k=0,n=16) 3.58e8, 6040.0, 2.0 e-28, 7.71e-5, 1.93e-30        # STAB U

# ONLY U FULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu_full...)
push!(wrc[2], ru_full...)
push!(wrc[3], cu_full...)
A = assemble_matrix(assem, wrc)
res_stab_ufull = compute_quantities(A,b,dΩ)
# (k=0,n=16) 3.58e8, 6050.0, 1.67e-27, 7.71e-5, 4.26e-30        # STAB U_FULL

# U + P + DIVU
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu...)
push!(wrc[2], ru...)
push!(wrc[3], cu...)
push!(wrc[1], wp...)
push!(wrc[2], rp...)
push!(wrc[3], cp...)
push!(wrc[1], wdiv...)
push!(wrc[2], rdiv...)
push!(wrc[3], cdiv...)
A = assemble_matrix(assem, wrc)
res_stab_updivu = compute_quantities(A,b,dΩ)
# (k=0,n=16) 7.59e6, 55500.0, 1.64e-26, 1.46e-4, 4.20e-29        # U + P + DIVU

# U_FULL + P + DIVU
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu_full...)
push!(wrc[2], ru_full...)
push!(wrc[3], cu_full...)
push!(wrc[1], wp...)
push!(wrc[2], rp...)
push!(wrc[3], cp...)
push!(wrc[1], wdiv...)
push!(wrc[2], rdiv...)
push!(wrc[3], cdiv...)
A = assemble_matrix(assem, wrc)
res_stab_ufullpdivu = compute_quantities(A,b,dΩ)
# (k=0,n=16) 7.59e6, 55500.0, 6.17e-26, 1.45e-4, 9.62e-29        # U_FULL + P + DIVU

# U_FULL + P + DIVU_FULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu_full...)
push!(wrc[2], ru_full...)
push!(wrc[3], cu_full...)
push!(wrc[1], wp...)
push!(wrc[2], rp...)
push!(wrc[3], cp...)
push!(wrc[1], wdiv_full...)
push!(wrc[2], rdiv_full...)
push!(wrc[3], cdiv_full...)
A = assemble_matrix(assem, wrc)
res_stab_ufullpdivufull = compute_quantities(A,b,dΩ)
# (k=0,n=16) 7.64e6, 54500.0, 8.79e-26, 1.46e-4, 1.98e-28       # U_FULL + P + DIVU_FULL

# U + P_FULL + DIVU
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu...)
push!(wrc[2], ru...)
push!(wrc[3], cu...)
push!(wrc[1], wp_full...)
push!(wrc[2], rp_full...)
push!(wrc[3], cp_full...)
push!(wrc[1], wdiv...)
push!(wrc[2], rdiv...)
push!(wrc[3], cdiv...)
A = assemble_matrix(assem, wrc)
res_stab_upfulldivu = compute_quantities(A,b,dΩ)
# (k=0,n=16) 1.28e7, 55500.0, 1.64e-26, 1.16e-4, 4.20e-29        # U + P_FULL + DIVU

# U_FULL + P_FULL + DIVU
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu_full...)
push!(wrc[2], ru_full...)
push!(wrc[3], cu_full...)
push!(wrc[1], wp_full...)
push!(wrc[2], rp_full...)
push!(wrc[3], cp_full...)
push!(wrc[1], wdiv...)
push!(wrc[2], rdiv...)
push!(wrc[3], cdiv...)
A = assemble_matrix(assem, wrc)
res_stab_ufullpfulldivu = compute_quantities(A,b,dΩ)
# (k=0,n=16) 1.28e7, 55500.0, 6.17e-26, 1.16e-4, 9.62e-29        # U_FULL + P_FULL + DIVU

# U_FULL + P_FULL + DIVU_FULL
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
A = assemble_matrix(assem, wrc)
res_stab_ufullpfulldivufull = compute_quantities(A,b,dΩ)
# (k=0,n=16) 1.29e7, 54500.0, 8.79e-26, 1.16e-4, 1.98e-28       # U_FULL + P_FULL + DIVU_FULL

# U + P + DIVU_FULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu...)
push!(wrc[2], ru...)
push!(wrc[3], cu...)
push!(wrc[1], wp...)
push!(wrc[2], rp...)
push!(wrc[3], cp...)
push!(wrc[1], wdiv_full...)
push!(wrc[2], rdiv_full...)
push!(wrc[3], cdiv_full...)
A = assemble_matrix(assem, wrc)
res_stab_updivufull = compute_quantities(A,b,dΩ)
# (k=0,n=16) 7.64e6, 54500.0, 2.75e-26, 1.46e-4, 1.14e-28       # U + P + DIVU_FULL



## Create collect function on cut cells only
"""
  dv, du: Test and trial basis functions. They may be components of a MultiFieldCellField

  # Compute and assemble the bulk penalty stabilization term
  # ∫( (div_dv-div_dv_l2_proj_agg_cells)*(div_du-div_du_l2_proj_agg_cells))*dΩ_cut_cells
"""
function div_penalty_stabilization_collect_cell_matrix_on_cut_cells_full(agg_cells_to_aggregate,
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

wdiv_full, rdiv_full, cdiv_full= div_penalty_stabilization_collect_cell_matrix_on_cut_cells_full(agg_cells_to_aggregate,
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
