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
nint = 2

# cut length
# ε    = 0.2   # not so small cut
ε    = 0.2e-2  # smaller cut
# ε    = 0.2e-6  # smallest cut

pmin = Point(0.0,0.0)
pmax = Point(1.0,1.0)
function setup_geometry(nint, ε, pmin, pmax)
   nbg     = nint + 2 + 2 # number of elements in the background mesh (2 from cut, 2 dummy)
   dp      = pmax - pmin
   h       = dp[1]/nint
   bgpmin  = pmin - Point(2*h,2*h)
   bgpmax  = pmax + Point(2*h,2*h)
   bgdp    = bgpmax - bgpmin
   partition = (nbg,nbg)
   bgmodel = CartesianDiscreteModel(bgpmin,bgpmax,partition)
   hbg     = bgdp[1]/nbg
   @assert abs(hbg - h) < 10e-10
   @assert abs(ε/h) < 0.5

   # The following is based on square (part of GridapEmbedded):
   e1 = VectorValue(1,0)
   e2 = VectorValue(0,1)
   x0 = pmin + Point(0.5*dp[1],0.5*dp[2])
   L1 = dp[1] + 2*ε
   L2 = dp[2] + 2*ε
   plane1 = plane(x0=x0-0.5*L2*e2,v=-e2,name="bottom")
   plane2 = plane(x0=x0+0.5*L1*e1,v= e1,name="right")
   plane3 = plane(x0=x0+0.5*L2*e2,v= e2,name="top")
   plane4 = plane(x0=x0-0.5*L1*e1,v=-e1,name="left")

   geo12 = intersect(plane1,plane2)
   geo34 = intersect(plane3,plane4)

   square = intersect(geo12,geo34)
   cutgeo = cut(bgmodel, square)

   # TODO: not necessary to export square as no cuts..
   bgmodel, cutgeo, h
end
bgmodel, cutgeo, h = setup_geometry(nint, ε, pmin, pmax)

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
degree=2*2*(order+1)
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
                                                        p_lhs,
                                                        P_Ωbg_cut_cell_dof_ids, 
                                                        P_agg_cells_local_dof_ids,
                                                        P_cut_cells_to_aggregate_dof_ids,
                                                        γ,
                                                        dΩbg_cut_cells)   

wdiv_full, rdiv_full, cdiv_full = div_penalty_stabilization_collect_cell_matrix_on_cut_cells_full(agg_cells_to_aggregate,
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

# MIXED STAB TERMS
wpmix, rpmix, cpmix = pmix_penalty_stabilization_collect_cell_matrix_on_cut_cells(agg_cells_to_aggregate,
                                                               ref_agg_cell_to_ref_bb_map,
                                                               dΩbg_agg_cells,
                                                               dq,    # Test basis
                                                               dp,    # Trial basis (to project)
                                                               qbb,  # Bounding box space test basis
                                                               p_lhs,
                                                               U_Ωbg_cut_cell_dof_ids,
                                                               P_Ωbg_cut_cell_dof_ids,
                                                               P_agg_cells_local_dof_ids,
                                                               P_cut_cells_to_aggregate_dof_ids,
                                                               γ,
                                                               dΩbg_cut_cells,
                                                               dv) 

wdmix, rdmix, cdmix = dmix_penalty_stabilization_collect_cell_matrix_on_cut_cells(agg_cells_to_aggregate,
                                                               ref_agg_cell_to_ref_bb_map,
                                                               dΩbg_agg_cells,
                                                               dq,    # Test basis
                                                               dp,    # Trial basis (to project)
                                                               qbb,  # Bounding box space test basis
                                                               p_lhs,
                                                               U_Ωbg_cut_cell_dof_ids,
                                                               P_Ωbg_cut_cell_dof_ids,
                                                               P_agg_cells_local_dof_ids,
                                                               P_cut_cells_to_aggregate_dof_ids,
                                                               γ,
                                                               dΩbg_cut_cells,
                                                               du) 

###########################################
###      VISUALIZATION                  ###
###########################################
writevtk(Ωbg,"trian-bg")
writevtk(Ωbg_cut_cells,"trian-bg-cut-cells")
colors = color_aggregates(aggregates,bgmodel)
writevtk(Ωbg,"trian-bg-with-coloured-aggregates", celldata=["aggregate"=>aggregates,"color"=>colors])
writevtk(Ωbg_agg_cells,"trian-bg-agg-cells")
Ωbg_int_cells=view(Ωbg,int_cells)
Ωbg_root_cells=view(Ωbg,root_cells)
writevtk(Ωbg_int_cells,"trian-bg-int-cells")
writevtk(Ωbg_root_cells,"trian-bg-root-cells")
writevtk(Ω,"trian-phys")
writevtk(Ωact,"trian-act")

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

function plot_quantities(A,b,Ω;filename="results")
   sol_x = A\b
   xh = FEFunction(X, sol_x)
   uh,ph = xh
   euh = uex-uh
   eph = pex-ph
   edivuh = divuex-(∇⋅uh)
   writevtk(Ω,filename,cellfields=["uh"=>uh,"ph"=>ph,"divuh"=>(∇⋅uh),"euh"=>euh,"eph"=>eph,"edivuh"=>edivuh])
end

## RHS STAB
# rhs_g = divuex # TODO: pass function on rather than constant value!!!!!
vecw_dmix, vecr_dmix = dmix_penalty_stabilization_collect_cell_vector_on_cut_cells(agg_cells_to_aggregate,
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
   dΩbg_cut_cells,
   divuex(1.0)) 

## TEST WITH MANUFACTURED SOLUTIONS
a((u,p),(v,q))=∫(v⋅u+q*p+(∇⋅u)*(∇⋅v))dΩ 
assem=SparseMatrixAssembler(X,Y)
l((v,q))=∫(v⋅uex+q*pex+(∇⋅v)*divuex)dΩ
b = assemble_vector(l, Y)

#RHS_DMIX
vec_wr=Gridap.FESpaces.collect_cell_vector(Y,l(dy))
push!(vec_wr[1],vecw_dmix...)
push!(vec_wr[2],vecr_dmix...)
vec_b = assemble_vector(assem, vec_wr)
# res_nostab = compute_quantities(A,vec_b,dΩ)

# NO STAB
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
A = assemble_matrix(assem, wrc)
res_nostab = compute_quantities(A,b,dΩ)

# ONLY U
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu...)
push!(wrc[2], ru...)
push!(wrc[3], cu...)
A = assemble_matrix(assem, wrc)
res_stab_u = compute_quantities(A,b,dΩ)

# ONLY U FULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu_full...)
push!(wrc[2], ru_full...)
push!(wrc[3], cu_full...)
A = assemble_matrix(assem, wrc)
res_stab_ufull = compute_quantities(A,b,dΩ)

# ONLY P 
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wp...)
push!(wrc[2], rp...)
push!(wrc[3], cp...)
A = assemble_matrix(assem, wrc)
res_stab_p = compute_quantities(A,b,dΩ)

# ONLY PFULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wp_full...)
push!(wrc[2], rp_full...)
push!(wrc[3], cp_full...)
A = assemble_matrix(assem, wrc)
res_stab_pfull = compute_quantities(A,b,dΩ)

# ONLY DIV
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wdiv...)
push!(wrc[2], rdiv...)
push!(wrc[3], cdiv...)
A = assemble_matrix(assem, wrc)
res_stab_div = compute_quantities(A,b,dΩ)

# ONLY DIVU FULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wdiv_full...)
push!(wrc[2], rdiv_full...)
push!(wrc[3], cdiv_full...)
A = assemble_matrix(assem, wrc)
res_stab_divfull = compute_quantities(A,b,dΩ)

# ONLY PMIX
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wpmix...)
push!(wrc[2], rpmix...)
push!(wrc[3], cpmix...)
A = assemble_matrix(assem, wrc)
# res_stab_pmix = compute_quantities(A,b,dΩ)
res_stab_pmix = compute_quantities(A,b,dΩ)

# ONLY DMIX
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wdmix...)
push!(wrc[2], rdmix...)
push!(wrc[3], cdmix...)
A = assemble_matrix(assem, wrc)
res_stab_dmix = compute_quantities(A,vec_b,dΩ)

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

# UFULL + PFULL + DIVUFULL
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

# U + DIVU + PMIX + DMIX
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu...)
push!(wrc[2], ru...)
push!(wrc[3], cu...)
push!(wrc[1], wdiv...)
push!(wrc[2], rdiv...)
push!(wrc[3], cdiv...)
push!(wrc[1], wpmix...)
push!(wrc[2], rpmix...)
push!(wrc[3], cpmix...)
push!(wrc[1], wdmix...)
push!(wrc[2], rdmix...)
push!(wrc[3], cdmix...)
A = assemble_matrix(assem, wrc)
res_stab_udivupmixdmix = compute_quantities(A,vec_b,dΩ)

# UFULL + DIVUFULL + PMIX + DMIX
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], wu_full...)
push!(wrc[2], ru_full...)
push!(wrc[3], cu_full...)
push!(wrc[1], wdiv_full...)
push!(wrc[2], rdiv_full...)
push!(wrc[3], cdiv_full...)
push!(wrc[1], wpmix...)
push!(wrc[2], rpmix...)
push!(wrc[3], cpmix...)
push!(wrc[1], wdmix...)
push!(wrc[2], rdmix...)
push!(wrc[3], cdmix...)
A = assemble_matrix(assem, wrc)
res_stab_ufulldivufullpmixdmix = compute_quantities(A,vec_b,dΩ)

## PRINT DATA:
# print(res_nostab, res_stab_u, res_stab_ufull, res_stab_p, res_stab_pfull, res_stab_div, res_stab_divfull, res_stab_pmix, res_stab_dmix, res_stab_updivu, res_stab_ufullpfulldivufull, res_stab_udivupmixdmix, res_stab_ufulldivufullpmixdmix)
res_nostab
res_stab_u
res_stab_ufull
res_stab_p
res_stab_pfull
res_stab_div
res_stab_divfull
res_stab_pmix
res_stab_dmix
res_stab_updivu
res_stab_ufullpfulldivufull
res_stab_udivupmixdmix
res_stab_ufulldivufullpmixdmix

# OLD:

############## (k=0,n=16) #####################################################################
# 2.31e9, 6040.0,  1.20e-27, 7.71e-5, 1.77e-30       # -- NO STAB --
# 3.58e8, 6040.0,  2.0 e-28, 7.71e-5, 1.93e-30       # U
# 3.58e8, 6050.0, 9.85e-28,  7.71e-5, 4.32e-30       # U_FULL
# 7.59e6, 55500.0, 1.14e-26, 1.46e-4, 1.38e-29       # U + P + DIVU
# 7.59e6, 55500.0, 6.17e-26, 1.45e-4, 9.62e-29       # U_FULL + P      + DIVU
# 7.64e6, 54500.0, 8.79e-26, 1.46e-4, 1.98e-28       # U_FULL + P      + DIVU_FULL (*)
# 1.28e7, 55500.0, 1.64e-26, 1.16e-4, 4.20e-29       # U      + P_FULL + DIVU
# 1.28e7, 55500.0, 6.17e-26, 1.16e-4, 9.62e-29       # U_FULL + P_FULL + DIVU
# 1.29e7, 54500.0, 8.79e-26, 1.16e-4, 1.98e-28       # U_FULL + P_FULL + DIVU_FULL
# 7.64e6, 54500.0, 2.75e-26, 1.46e-4, 1.14e-28       # U      + P      + DIVU_FULL (*)
###############################################################################################

############## (k=2,n=16) #####################################################################
# 2.82e29, 9.21e8, 4.0e-12, 7.79e-28, 2.34e-13        # -- NO STAB --
# 2.82e29, 9.21e8, 2.05e-18, 7.79e-28, 1.66e-19       # U
# 2.82e29, 9.21e8, 1.23e-16, 7.79e-28, 1.39e-18       # U_FULL
# 4.72e16, 9.45e9, 3.08e-14, 4.07e-29, 1.76e-17       # U      + P      + DIVU
# 4.72e16, 9.45e9, 2.52e-11, 4.07e-29, 3.17e-16       # U_FULL + P      + DIVU
# 3.00e18, 9.45e9, 3.08e-14, 1.65e-28, 1.76e-17       # U      + P_FULL + DIVU
# 3.00e18, 9.45e9, 2.52e-11, 1.65e-28, 3.17e-16       # U_FULL + P_FULL + DIVU
###############################################################################################