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
                                                        
## TEST WITH MANUFACTURED SOLUTIONS
a((u,p),(v,q))=∫(v⋅u+q*p+(∇⋅u)*(∇⋅v))dΩ 
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
assem=SparseMatrixAssembler(X,Y)
A_nostab = assemble_matrix(assem, wrc)
l((v,q))=∫(v⋅uex+q*pex+(∇⋅v)*divuex)dΩ
b = assemble_vector(l, Y)

push!(wrc[1], wu...)
push!(wrc[2], ru...)
push!(wrc[3], cu...)

push!(wrc[1], wp...)
push!(wrc[2], rp...)
push!(wrc[3], cp...)

A_stabup = assemble_matrix(assem, wrc) 
cond_nostab = cond(Array(A_nostab))
cond_stabup = cond(Array(A_stabup))

# CHECK errors when only the p and u stabilization terms are included.
global_l2_proj_dofs_upstab = A_stabup\b
xh_upstab = FEFunction(X, global_l2_proj_dofs_upstab)
uh_upstab,ph_upstab = xh_upstab
euh_upstab = uex-uh_upstab
eph_upstab = pex-ph_upstab
edivuh_upstab = divuex-(∇⋅uh_upstab)
norm_euh_upstab = sum(∫(euh_upstab⋅euh_upstab)*dΩ)
norm_eph_upstab = sum(∫(eph_upstab*eph_upstab)*dΩ)
norm_edivuh_upstab = sum(∫(edivuh_upstab⋅edivuh_upstab)*dΩ)
# No stabilization (n=16) (k=0) κ = 2.3e9 vs (k=2) κ = 2.8e29
# k=0 (n=16) yields euh = 2.0e-28, eph = 0.00015, edivuh = 1.9e-30, κ = 9.2e5
# k=2 (n=16) yields euh = 2.0e-18, eph = 4.0e-29, edivuh = 1.7e-19, κ = 4.4e15

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
push!(wrc[1], wdiv...)
push!(wrc[2], rdiv...)
push!(wrc[3], cdiv...)

A_stab = assemble_matrix(assem, wrc) 
cond_stab = cond(Array(A_stab))

global_l2_proj_dofs_stab = A_stab\b
xh_stab = FEFunction(X, global_l2_proj_dofs_stab)
uh_stab,ph_stab = xh_stab
euh_stab = uex-uh_stab
eph_stab = pex-ph_stab
edivuh_stab = divuex-(∇⋅uh_stab)
norm_euh_stab = sum(∫(euh_stab⋅euh_stab)*dΩ)
norm_eph_stab = sum(∫(eph_stab*eph_stab)*dΩ)
norm_edivuh_stab = sum(∫(edivuh_stab⋅edivuh_stab)*dΩ)
# k=0 (n=16) yields euh = 1.6e-26, eph = 0.00015, edivuh = 4.2e-29, κ = 7.6e6
# k=2 (n=16) yields euh = 3.1e-14, eph = 4.1e-29, edivuh = 1.8e-17, κ = 4.7e16