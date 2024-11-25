using Gridap
using GridapEmbedded
using FillArrays
using LinearAlgebra

include("aggregates_bounding_boxes_tools.jl")
include("bulk_ghost_penalty_stab_tools.jl")
include("fields_and_blocks_tools.jl")
include("BulkGhostPenaltyAssembleMaps.jl")

# Manufactured solution
order = 2
uex(x) = -VectorValue(2*x[1],2*x[2])
pex(x) = (x[1]^2 + x[2]^2)
divuex(x) = -4.0
# Select geometry
R = 0.2
geom = disk(R, x0=Point(0.5,0.5))

# Setup background model
n=10
partition = (n,n)
box = get_metadata(geom)
bgmodel = CartesianDiscreteModel((0,1,0,1),partition)
dp = box.pmax - box.pmin
h = dp[1]/n

# Cut the background model with the mesh
cutdisk = cut(bgmodel,geom)

# Compute mapping among background model 
# cut cells and interior cells 
strategy = AggregateAllCutCells()
aggregates = aggregate(strategy,cutdisk)


aggregate_to_cells=setup_aggregate_to_cells(aggregates)
aggregates_bounding_box_model=
       setup_aggregates_bounding_box_model(bgmodel,aggregate_to_cells)

# colors = color_aggregates(aggregates,bgmodel)
# writevtk(Triangulation(bgmodel),"trian",celldata=["cellin"=>aggregates,"color"=>colors])       
# writevtk(aggregates_bounding_box_model, "bb_model")
# writevtk(bgmodel, "bg_model")

# Set up objects required to compute both LHS and RHS of the L2 projection

# Set up global spaces 
Ωhact = Triangulation(cutdisk,ACTIVE)

V = FESpace(Ωhact, ReferenceFE(raviart_thomas,Float64,order),conformity=:HDiv)
Q = FESpace(Ωhact, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
U = TrialFESpace(V)
P = TrialFESpace(Q)
Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

agg_cells=setup_agg_cells(aggregate_to_cells)

# ref_agg_cell_to_agg_cell_map: \hat{K} -> K 
Ω=Triangulation(bgmodel)
Ωagg_cells=view(Ω,agg_cells)
ref_agg_cell_to_agg_cell_map=get_cell_map(Ωagg_cells)

agg_cells_to_aggregate=setup_agg_cells_to_aggregate(aggregate_to_cells)

ref_agg_cell_to_ref_bb_map=setup_ref_agg_cell_to_ref_bb_map(aggregates_bounding_box_model,
                                                            agg_cells_to_aggregate)

# Compute LHS of L2 projection 
degree=2*(order+1)
dΩagg_cells = Measure(Ωagg_cells,degree)
reffe =ReferenceFE(lagrangian,Float64,order) # Here we MUST use a Q space (not a P space!)
Qbb=FESpace(aggregates_bounding_box_model,reffe,conformity=:L2) # We need a DG space to represent the L2 projection
Pbb=TrialFESpace(Qbb)
pbb=get_trial_fe_basis(Pbb)
qbb=get_fe_basis(Qbb)

reffe=ReferenceFE(raviart_thomas,Float64,order)
Vbb=FESpace(aggregates_bounding_box_model,reffe,conformity=:L2)
Ubb=TrialFESpace(Vbb)
ubb=get_trial_fe_basis(Ubb)
vbb=get_fe_basis(Vbb)

aggregate_to_local_cells=setup_aggregate_to_local_cells(aggregate_to_cells)

p_lhs=set_up_bulk_ghost_penalty_lhs(aggregates_bounding_box_model,
                                  agg_cells_to_aggregate,
                                  ref_agg_cell_to_ref_bb_map,
                                  dΩagg_cells,
                                  pbb,
                                  qbb)

u_lhs=set_up_bulk_ghost_penalty_lhs(aggregates_bounding_box_model,
                                  agg_cells_to_aggregate,
                                  ref_agg_cell_to_ref_bb_map,
                                  dΩagg_cells,
                                  ubb,
                                  vbb)                             

# Compute contributions to the RHS of the L2 projection
dx    = get_trial_fe_basis(X)
dy    = get_fe_basis(Y)
du,dp = dx  
dv,dq = dy 

Ωagg_cell_dof_ids   = get_cell_dof_ids(X,Ωagg_cells)
U_Ωagg_cell_dof_ids = _restrict_to_block(Ωagg_cell_dof_ids, 1) 
P_Ωagg_cell_dof_ids = _restrict_to_block(Ωagg_cell_dof_ids, 2)
# Ωagg_cell_dof_ids[1][1] # cell 1, velo dof
# Ωagg_cell_dof_ids[1][2] # cell 1, pressure dofs
 
### BEGIN TESTING CODE
# This code is just for testing purposes, so I have commented it out 
# It allows to evaluate the LHS of the L2 projection corresponding to a 
# particular FE function, instead of a basis 
# uhex=interpolate(uex,Ustd)
# agg_cells_lhs_contribs_uhex=get_array(∫(vbb_Ωagg_cells*uhex)dΩagg_cells)
# ass_lhs_map_uhex=AssembleLhsMap(agg_cells_lhs_contribs_uhex)
# lhs_uhex=lazy_map(ass_lhs_map_uhex,aggregate_to_local_cells)
### END TESTING CODE


U_agg_cells_local_dof_ids=
   compute_agg_cells_local_dof_ids(U_Ωagg_cell_dof_ids, aggregate_to_local_cells)
P_agg_cells_local_dof_ids=
   compute_agg_cells_local_dof_ids(P_Ωagg_cell_dof_ids, aggregate_to_local_cells)


function set_up_h_U(aggregates_bounding_box_model,
                    agg_cells_to_aggregate,
                    Ωagg_cells)
    degree = 0 # We are integrating a constant function 
               # Thus, degree=0 is enough for exact integration
    Ωbb  = Triangulation(aggregates_bounding_box_model)
    dΩbb = Measure(Ωbb, degree)
    h_U_array = get_array(∫(1.0)dΩbb)
    h_U_array = lazy_map(Reindex(h_U_array), agg_cells_to_aggregate)
    CellField(h_U_array, Ωagg_cells)
end 

U_aggregate_dof_ids=compute_aggregate_dof_ids(U_Ωagg_cell_dof_ids,aggregate_to_cells)
U_agg_cells_to_aggregate_dof_ids=lazy_map(Reindex(U_aggregate_dof_ids),agg_cells_to_aggregate)

P_aggregate_dof_ids=compute_aggregate_dof_ids(P_Ωagg_cell_dof_ids,aggregate_to_cells)
P_agg_cells_to_aggregate_dof_ids=lazy_map(Reindex(P_aggregate_dof_ids),agg_cells_to_aggregate)

γ = 10.0 # Interior bulk-penalty stabilization parameter
         # (@amartinhuertas no idea what a reasonable value is)

h_U = set_up_h_U(aggregates_bounding_box_model, agg_cells_to_aggregate, Ωagg_cells)         

# Manually set up the arrays that collect_cell_matrix would return automatically
wp,rp,cp=bulk_ghost_penalty_stabilization_collect_cell_matrix(agg_cells_to_aggregate,
                                                        ref_agg_cell_to_ref_bb_map,
                                                        dΩagg_cells,
                                                        dq,    # Test basis
                                                        dp,    # Trial basis (to project)
                                                        qbb,  # Bounding box space test basis
                                                        p_lhs,
                                                        P_Ωagg_cell_dof_ids, 
                                                        P_agg_cells_local_dof_ids,
                                                        P_agg_cells_to_aggregate_dof_ids,
                                                        h_U,
                                                        γ)

wu,ru,cu=bulk_ghost_penalty_stabilization_collect_cell_matrix(agg_cells_to_aggregate,
                                                        ref_agg_cell_to_ref_bb_map,
                                                        dΩagg_cells,
                                                        dv,    # Test basis
                                                        du,    # Trial basis (to project)
                                                        vbb,  # Bounding box space test basis
                                                        u_lhs,
                                                        U_Ωagg_cell_dof_ids, 
                                                        U_agg_cells_local_dof_ids,
                                                        U_agg_cells_to_aggregate_dof_ids,
                                                        h_U,
                                                        γ)

wdivu,rdivu,cdivu=div_penalty_stabilization_collect_cell_matrix(agg_cells_to_aggregate,
                                                        ref_agg_cell_to_ref_bb_map,
                                                        dΩagg_cells,
                                                        dv,    # Test basis
                                                        du,    # Trial basis (to project)
                                                        qbb,  # Bounding box space test basis
                                                        p_lhs,
                                                        U_Ωagg_cell_dof_ids, 
                                                        U_agg_cells_local_dof_ids,
                                                        U_agg_cells_to_aggregate_dof_ids,
                                                        h_U,
                                                        γ)

Ωcut = Triangulation(cutdisk,PHYSICAL)
dΩcut = Measure(Ωcut,degree)


## TEST 1
a((u,p),(v,q))=∫(v⋅u+q*p+(∇⋅u)*(∇⋅v))dΩcut 
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))

assem=SparseMatrixAssembler(X,Y)
Anostab=assemble_matrix(assem, wrc)
cond_nostab = cond(Array(Anostab))

# Set up rhs global projection 
l((v,q))=∫(v⋅uex+q*pex+(∇⋅v)*divuex)dΩcut
b = assemble_vector(l, Y)

## U + P STABILIZATION 
# Add the bulk penalty stabilization terms for u and p to wrc
push!(wrc[1], wp...)
push!(wrc[2], rp...)
push!(wrc[3], cp...)

push!(wrc[1], wu...)
push!(wrc[2], ru...)
push!(wrc[3], cu...)

# Check the condition number when only the p and u stabilization terms are included.
Aupstab=assemble_matrix(assem, wrc)
cond_upstab = cond(Array(Aupstab))

# CHECK errors when only the p and u stabilization terms are included.
global_l2_proj_dofs_upstab = Aupstab\b
xh_upstab = FEFunction(X, global_l2_proj_dofs_upstab)
uh_upstab,ph_upstab = xh_upstab
euh_upstab = uex-uh_upstab
eph_upstab = pex-ph_upstab
edivuh_upstab = divuex-(∇⋅uh_upstab)
norm_euh_upstab = sum(∫(euh_upstab⋅euh_upstab)*dΩcut)
norm_eph_upstab = sum(∫(eph_upstab*eph_upstab)*dΩcut)
norm_edivuh_upstab = sum(∫(edivuh_upstab⋅edivuh_upstab)*dΩcut)
# No stabilization (k=0) κ = 1.1e35 vs (k=2) κ = Inf
# k=0 (n=10) yields euh = 3.0e-28, eph = 0.00055, edivuh = 1.3e-30, κ = 1.4e5
# k=2 (n=10) yields euh = 1.8e-18, eph = 1.7e-28, edivuh = 7.4e-21, κ= 1.7e13

## Now add div stabilization to A
push!(wrc[1], wdivu...)
push!(wrc[2], rdivu...)
push!(wrc[3], cdivu...)

# Assembly 
Awithstab=assemble_matrix(assem, wrc)
cond_withstab = cond(Array(Awithstab))

global_l2_proj_dofs_withstab = Awithstab\b
xh_withstab = FEFunction(X, global_l2_proj_dofs_withstab)
uh_withstab,ph_withstab = xh_withstab
euh_withstab = uex-uh_withstab
eph_withstab = pex-ph_withstab
edivuh_withstab = divuex-(∇⋅uh_withstab)
norm_euh_withstab = sum(∫(euh_withstab⋅euh_withstab)*dΩcut)
norm_eph_withstab = sum(∫(eph_withstab*eph_withstab)*dΩcut)
norm_edivuh_withstab = sum(∫(edivuh_withstab⋅edivuh_withstab)*dΩcut)
# k=0 (n=10) yields euh = 1.2e-26, eph = 0.00055, edivuh = 7.0e-29, κ = 1.7e6
# k=2 (n=10) yields euh = 1.6e-15, eph = 1.7e-28, edivuh = 7.1e-19, κ= 1.8e14