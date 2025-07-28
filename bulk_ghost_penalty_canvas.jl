using Gridap
using GridapEmbedded

include("./src/BGP/BGP.jl")

# Problem selection
problem = 1 # 0 = Manufactured solution (L2-like projection), 1 = Darcy problem 

# Manufactured solution
order  = 0
uex(x) = -VectorValue(2*x[1],2*x[2])
pex(x) = (x[1]^2 + x[2]^2)
divuex(x) = -4.0

# Select geometry
nint = 3    # number of (uncut) interior elements along single direction 
# cut length
ε    = 0.2/2.0     # not so small cut
# ε    = 0.2e-2/2.0  # smaller cut
# ε    = 0.2e-6/2.0  # smallest cut
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

   bgmodel, cutgeo, h
end
bgmodel, cutgeo, h= setup_geometry(nint, ε, pmin, pmax)

# Setup aggregates
strategy  = AggregateAllCutCells()
aggregates= aggregate(strategy,cutgeo)
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

# Setup agg cells and triangulation
agg_cells    =flatten(aggregate_to_cells) #[CLEAN]: replaced setup_agg_cells
Ωbg_agg_cells=view(Ωbg,agg_cells)

# Pressure space
if problem==0
   Q = FESpace(Ωact, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
elseif problem==1
   # Set up zero-mean pressure space, with fixed interior dof
   Qnzm = FESpace(Ωact, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
   int_cells = restrict_cells(cutgeo,IN) # global (background mesh') cell identifiers of interior cells 
   nonagg_int_cell = int_cells[findfirst(!in(agg_cells),int_cells)]   # cell identifier (background mesh) of first interior cell not in aggregate
   local_id = findfirst(isequal(nonagg_int_cell),Ωact.tface_to_mface) # local (active mesh') cell id of interior cell not in aggregate
   dof_to_fix = get_cell_dof_ids(Qnzm)[local_id][1]
   spaceWithConstantFixed = Gridap.FESpaces.FESpaceWithConstantFixed(Qnzm,true,Int64(dof_to_fix)) 
   Qzm_vol_i = assemble_vector(v->∫(v)*dΩ,Qnzm)
   Qzm_vol = sum(Qzm_vol_i)
   Q     = Gridap.FESpaces.ZeroMeanFESpace(spaceWithConstantFixed,Qzm_vol_i,Qzm_vol)
end

# Set up global spaces 
V = FESpace(Ωact, ReferenceFE(raviart_thomas,Float64,order),conformity=:HDiv)
U = TrialFESpace(V)
P = TrialFESpace(Q)
Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])
dx    = get_trial_fe_basis(X)
dy    = get_fe_basis(Y)
du,dp = dx  
dv,dq = dy 

# ref_agg_cell_to_agg_cell_map: \hat{K} -> K 
ref_agg_cell_to_agg_cell_map=get_cell_map(Ωbg_agg_cells)
agg_cells_to_aggregate      =setup_cells_to_aggregate(aggregate_to_cells)
ref_agg_cell_to_ref_bb_map  =setup_ref_agg_cell_to_ref_bb_map(aggregates_bounding_box_model,
                                                            agg_cells_to_aggregate,ref_agg_cell_to_agg_cell_map)

# Spaces on bounding boxes                                                            
reffeₚ_bb =ReferenceFE(lagrangian,Float64,order) 
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

# Selecting relevant global dofs ids of aggregate cells (from background mesh)
Ωbg_agg_cell_dof_ids   = get_cell_dof_ids(X,Ωbg_agg_cells)
U_Ωbg_agg_cell_dof_ids = _restrict_to_block(Ωbg_agg_cell_dof_ids, 1) 
P_Ωbg_agg_cell_dof_ids = _restrict_to_block(Ωbg_agg_cell_dof_ids, 2)

# Computing local (per aggregate) dof ids 
aggregate_to_local_cells=setup_aggregate_to_local_cells(aggregate_to_cells)
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

###########################################
### STABILIZATION ON Ωagg\Troot ###
###########################################

# Setup cut cells and triangulation
aggregate_to_cut_cells = restrict_aggregate_to_cells(cutgeo,aggregate_to_cells,GridapEmbedded.Interfaces.CUT)
cut_cells = flatten(aggregate_to_cut_cells)
#TO-DO: look into why cut_cells  = restrict_cells(cutgeo,GridapEmbedded.Interfaces.CUT) can not be used instead of the above
Ωbg_cut_cells   = view(Ωbg,cut_cells)
dΩbg_cut_cells  = Measure(Ωbg_cut_cells,degree)

# Selecting relevant global dofs ids of cut cells (from background mesh)
Ωbg_cut_cell_dof_ids   = get_cell_dof_ids(X,Ωbg_cut_cells)
U_Ωbg_cut_cell_dof_ids = _restrict_to_block(Ωbg_cut_cell_dof_ids, 1) 
P_Ωbg_cut_cell_dof_ids = _restrict_to_block(Ωbg_cut_cell_dof_ids, 2)

# Compute global dofs ids per aggregate and reindex these 
cut_cells_to_aggregate = setup_cells_to_aggregate(aggregate_to_cut_cells)
U_cut_cells_to_aggregate_dof_ids=lazy_map(Reindex(U_aggregate_dof_ids),cut_cells_to_aggregate)
P_cut_cells_to_aggregate_dof_ids=lazy_map(Reindex(P_aggregate_dof_ids),cut_cells_to_aggregate)

###########################################
### Setup projections
###########################################

du_proj_Vbb, dv_proj_Vbb = setup_L2_proj_in_bb_space(
   dΩbg_agg_cells,             # measure of aggregated cells in background domain
   ref_agg_cell_to_ref_bb_map, # map
   agg_cells_to_aggregate,     # 
   aggregate_to_local_cells,   # 
   du,                         # Trial basis (to project) 
   dv,                         # Test basis 
   ubb,                        # Trial basis of bounding box space Vbb
   vbb,                        # Test basis of bounding box space Vbb
   identity,                   # operation to be applied to u and v
   U_agg_cells_local_dof_ids)  # aggregates local dof ids for space U   

dp_proj_Qbb, dq_proj_Qbb = setup_L2_proj_in_bb_space(
   dΩbg_agg_cells,             # measure of aggregated cells in background domain
   ref_agg_cell_to_ref_bb_map, # map
   agg_cells_to_aggregate,     # 
   aggregate_to_local_cells,   # 
   dp,                         # Trial basis (to project) 
   dq,                         # Test basis 
   pbb,                        # Trial basis of bounding box space Qbb
   qbb,                        # Test basis of bounding box space Qbb
   identity,                   # operation to be applied to u and v
   P_agg_cells_local_dof_ids)  # aggregates local dof ids for space P   

div_du_proj_Qbb, div_dv_proj_Qbb = setup_L2_proj_in_bb_space(
   dΩbg_agg_cells,             # measure of aggregated cells in background domain
   ref_agg_cell_to_ref_bb_map, # map
   agg_cells_to_aggregate,     # 
   aggregate_to_local_cells,   # 
   du,                         # Trial basis (to project) 
   dv,                         # Test basis 
   pbb,                        # Trial basis of bounding box space Vbb
   qbb,                        # Test basis of bounding box space Vbb
   divergence,                 # operation to be applied to u and v
   U_agg_cells_local_dof_ids)  # aggregates local dof ids for space U   

# ###########################################
# ### BlockMap PREP STEPS #TODO MOVE INTO collect_cell_matrix_on_D?
# ###########################################
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

##########################################
### Setup stabilization terms (DIFF)   ###
##########################################
w_u_diff, r_u_diff, c_u_diff = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dΩbg_cut_cells,
   Ωbg_agg_cells,
   γ,
   du,
   dv,
   du_proj_Vbb,
   dv_proj_Vbb,
   U_Ωbg_cut_cell_dof_ids,
   U_Ωbg_cut_cell_dof_ids,
   U_cut_cells_to_aggregate_dof_ids,
   U_cut_cells_to_aggregate_dof_ids,
   identity,
   identity)

w_p_diff, r_p_diff, c_p_diff = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dΩbg_cut_cells,
   Ωbg_agg_cells,
   γ,
   dp,
   dq,
   dp_proj_Qbb,
   dq_proj_Qbb,
   P_Ωbg_cut_cell_dof_ids,
   P_Ωbg_cut_cell_dof_ids,
   P_cut_cells_to_aggregate_dof_ids,
   P_cut_cells_to_aggregate_dof_ids,
   identity,
   identity)

w_divu_diff, r_divu_diff, c_divu_diff = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dΩbg_cut_cells,
   Ωbg_agg_cells,
   γ,
   du,
   dv,
   div_du_proj_Qbb,
   div_dv_proj_Qbb,
   U_Ωbg_cut_cell_dof_ids,
   U_Ωbg_cut_cell_dof_ids,
   U_cut_cells_to_aggregate_dof_ids,
   U_cut_cells_to_aggregate_dof_ids,
   divergence,
   divergence)

w_divuq_diff, r_divuq_diff, c_divuq_diff = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dΩbg_cut_cells,
   Ωbg_agg_cells,
   γ,
   du,
   dq,
   div_du_proj_Qbb,
   dq_proj_Qbb,
   U_Ωbg_cut_cell_dof_ids,
   P_Ωbg_cut_cell_dof_ids,
   U_cut_cells_to_aggregate_dof_ids,
   P_cut_cells_to_aggregate_dof_ids,
   divergence,
   identity)

w_pdivv_diff, r_pdivv_diff, c_pdivv_diff = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dΩbg_cut_cells,
   Ωbg_agg_cells,
   γ,
   dp,
   dv,
   dp_proj_Qbb,
   div_dv_proj_Qbb,
   P_Ωbg_cut_cell_dof_ids,
   U_Ωbg_cut_cell_dof_ids,
   P_cut_cells_to_aggregate_dof_ids,
   U_cut_cells_to_aggregate_dof_ids,
   identity,
   divergence)

##########################################
### Setup stabilization terms (FULL)   ###
##########################################
w_u_full, r_u_full, c_u_full = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dΩbg_cut_cells,
   Ωbg_agg_cells,
   γ,
   du,
   dv,
   du_proj_Vbb,
   U_Ωbg_cut_cell_dof_ids,
   U_Ωbg_cut_cell_dof_ids,
   U_cut_cells_to_aggregate_dof_ids,
   identity,
   identity)

w_p_full, r_p_full, c_p_full = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dΩbg_cut_cells,
   Ωbg_agg_cells,
   γ,
   dp,
   dq,
   dp_proj_Qbb,
   P_Ωbg_cut_cell_dof_ids,
   P_Ωbg_cut_cell_dof_ids,
   P_cut_cells_to_aggregate_dof_ids,
   identity,
   identity)

w_divu_full, r_divu_full, c_divu_full = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dΩbg_cut_cells,
   Ωbg_agg_cells,
   γ,
   du,
   dv,
   div_du_proj_Qbb,
   U_Ωbg_cut_cell_dof_ids,
   U_Ωbg_cut_cell_dof_ids,
   U_cut_cells_to_aggregate_dof_ids,
   divergence,
   divergence)

w_divuq_full, r_divuq_full, c_divuq_full = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dΩbg_cut_cells,
   Ωbg_agg_cells,
   γ,
   du,
   dq,
   div_du_proj_Qbb,
   U_Ωbg_cut_cell_dof_ids,
   P_Ωbg_cut_cell_dof_ids,
   U_cut_cells_to_aggregate_dof_ids,
   divergence,
   identity)

w_pdivv_full, r_pdivv_full, c_pdivv_full = bulk_ghost_penalty_stabilization_collect_cell_matrix_on_D(dΩbg_cut_cells,
   Ωbg_agg_cells,
   γ,
   dp,
   dv,
   dp_proj_Qbb,
   P_Ωbg_cut_cell_dof_ids,
   U_Ωbg_cut_cell_dof_ids,
   P_cut_cells_to_aggregate_dof_ids,
   identity,
   divergence)

##########################################
### Setup rhs stabilization            ###
##########################################
rhs_func = divuex
rhs_func_proj_Qbb = setup_L2_proj_in_bb_space(dΩbg_agg_cells,
   ref_agg_cell_to_ref_bb_map, 
   agg_cells_to_aggregate,      
   aggregate_to_local_cells,    
   rhs_func,                   
   pbb,                        
   qbb) 

w_rhsq_diff, r_rhsq_diff = bulk_ghost_penalty_stabilization_collect_cell_vector_on_D(dΩbg_cut_cells,
   γ,  
   dq,
   dq_proj_Qbb,
   P_Ωbg_cut_cell_dof_ids,
   P_cut_cells_to_aggregate_dof_ids,
   identity,
   rhs_func,
   rhs_func_proj_Qbb)

w_rhsq_full, r_rhsq_full = bulk_ghost_penalty_stabilization_collect_cell_vector_on_D(dΩbg_cut_cells,
   γ,  
   dq, 
   P_Ωbg_cut_cell_dof_ids,
   identity,
   rhs_func,
   rhs_func_proj_Qbb)

## WEAK FORM
if problem==0
   a((u,p),(v,q))=∫(v⋅u+q*p+(∇⋅u)*(∇⋅v))dΩ 
   l((v,q))=∫(v⋅uex+q*pex+(∇⋅v)*divuex)dΩ
elseif problem==1
   ΓN    = EmbeddedBoundary(cutgeo)
   dΓN   = Measure(ΓN,degree)
   nN    = get_normal_vector(ΓN)
   s     = -1 
   β₀    = 1e3
   β     = β₀*h.^(s)
   m     =1
   uNeu = uex
   g    = divuex
   a((u,p), (v,q)) = ∫(v⋅u - (∇⋅v)*p - (∇⋅u)*q)dΩ + ∫(β*(u⋅nN)*(v⋅nN) + (v⋅nN)*p + m*(u⋅nN)*q)dΓN    
   l((v,q))        = ∫(- q*g)dΩ + ∫(β*(v⋅nN)*(uNeu⋅nN) + m*q*(uNeu⋅nN))dΓN 
else
   print("Problem not implemented. Try again with problem=0 (L2-like projection test) or problem=1 (Darcy test)")
end
assem=SparseMatrixAssembler(X,Y)
b = assemble_vector(l, Y)

# RHS
vec_rhsq_diff=Gridap.FESpaces.collect_cell_vector(Y,l(dy))
push!(vec_rhsq_diff[1],w_rhsq_diff...)
push!(vec_rhsq_diff[2],r_rhsq_diff...)
b_rhsq_diff = assemble_vector(assem, vec_rhsq_diff)

vec_rhsq_full=Gridap.FESpaces.collect_cell_vector(Y,l(dy))
push!(vec_rhsq_full[1],w_rhsq_full...)
push!(vec_rhsq_full[2],r_rhsq_full...)
b_rhsq_full = assemble_vector(assem, vec_rhsq_full)

function compute_quantities(problem,A,b,dΩ)
   cond_A =  cond(Array(A))
   norm_A =  norm(A)
   sol_x = A\b
   xh = FEFunction(X, sol_x)
   uh,ph = xh
   if problem==1
      area = sum(∫(1.0)dΩ)
      mean_p  = sum(∫(pex)dΩ)/area # mean presure exact sol
      ph = ph + mean_p
   end
   euh = uex-uh
   eph = pex-ph
   edivuh = divuex-(∇⋅uh)
   norm_euh = sum(∫(euh⋅euh)*dΩ)
   norm_eph = sum(∫(eph*eph)*dΩ)
   norm_edivuh = sum(∫(edivuh⋅edivuh)*dΩ)
   return round(cond_A,sigdigits=3), round(norm_A,sigdigits=3), round(norm_euh,sigdigits=3), round(norm_eph,sigdigits=3), round(norm_edivuh,sigdigits=3)
end

####### SOLVING TEST PROBLEMS (DIFF)

# NO STAB DIFF
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
A = assemble_matrix(assem, wrc)
res_nostab = compute_quantities(problem,A,b,dΩ)

# ONLY U DIFF
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_u_diff...)
push!(wrc[2], r_u_diff...)
push!(wrc[3], c_u_diff...)
A = assemble_matrix(assem, wrc)
res_stab_u_diff = compute_quantities(problem,A,b,dΩ)

# ONLY P DIFF
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_p_diff...)
push!(wrc[2], r_p_diff...)
push!(wrc[3], c_p_diff...)
A = assemble_matrix(assem, wrc)
res_stab_p_diff = compute_quantities(problem,A,b,dΩ)

# ONLY DIVU DIFF
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_divu_diff...)
push!(wrc[2], r_divu_diff...)
push!(wrc[3], c_divu_diff...)
A = assemble_matrix(assem, wrc)
res_stab_divu_diff = compute_quantities(problem,A,b,dΩ)

# ONLY DIVUQ DIFF
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_divuq_diff...)
push!(wrc[2], r_divuq_diff...)
push!(wrc[3], c_divuq_diff...)
A = assemble_matrix(assem, wrc)
res_stab_divuq_diff = compute_quantities(problem,A,b,dΩ)
res_stab_divuq_diff_rhs = compute_quantities(problem,A,b_rhsq_diff,dΩ)

# ONLY PDIVV DIFF
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_pdivv_diff...)
push!(wrc[2], r_pdivv_diff...)
push!(wrc[3], c_pdivv_diff...)
A = assemble_matrix(assem, wrc)
res_stab_pdivv_diff = compute_quantities(problem,A,b,dΩ)

####### SOLVING TEST PROBLEMS (FULL)

# ONLY U FULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_u_full...)
push!(wrc[2], r_u_full...)
push!(wrc[3], c_u_full...)
A = assemble_matrix(assem, wrc)
res_stab_u_full = compute_quantities(problem,A,b,dΩ)

# ONLY P FULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_p_full...)
push!(wrc[2], r_p_full...)
push!(wrc[3], c_p_full...)
A = assemble_matrix(assem, wrc)
res_stab_p_full = compute_quantities(problem,A,b,dΩ)

# ONLY DIVU FULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_divu_full...)
push!(wrc[2], r_divu_full...)
push!(wrc[3], c_divu_full...)
A = assemble_matrix(assem, wrc)
res_stab_divu_full = compute_quantities(problem,A,b,dΩ)

# ONLY DIVUQ FULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_divuq_full...)
push!(wrc[2], r_divuq_full...)
push!(wrc[3], c_divuq_full...)
A = assemble_matrix(assem, wrc)
res_stab_divuq_full = compute_quantities(problem,A,b,dΩ)
res_stab_divuq_full_rhs = compute_quantities(problem,A,b_rhsq_full,dΩ)

# ONLY PDIVV FULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_pdivv_full...)
push!(wrc[2], r_pdivv_full...)
push!(wrc[3], c_pdivv_full...)
A = assemble_matrix(assem, wrc)
res_stab_pdivv_full = compute_quantities(problem,A,b,dΩ)

## GENERATE EXCEL SHEET RESULTS``
print("results")
res_nostab
res_stab_u_full
res_stab_u_diff
res_stab_p_full
res_stab_p_diff
res_stab_divu_full
res_stab_divu_diff
res_stab_pdivv_full # THIS IS OK. 
res_stab_divuq_full_rhs # ERRORS SEEM TO BE OF SAME ORDER NOW.

# Alternatives for res_stab_divuq_full_rhs
res_stab_divuq_full #THIS SEEMS BETTER
res_stab_divuq_diff #not too relevant 
res_stab_divuq_full_rhs # seems incorrect
res_stab_divuq_diff_rhs # seems incorrect

# Lastly compare these:
res_stab_pdivv_full 
res_stab_pdivv_diff 

###
# UFULL + PFULL + DIVUFULL
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_u_full...)
push!(wrc[2], r_u_full...)
push!(wrc[3], c_u_full...)
push!(wrc[1], w_p_full...)
push!(wrc[2], r_p_full...)
push!(wrc[3], c_p_full...)
push!(wrc[1], w_divu_full...)
push!(wrc[2], r_divu_full...)
push!(wrc[3], c_divu_full...)
A = assemble_matrix(assem, wrc)
res_stab_updivu_full = compute_quantities(problem,A,b,dΩ) #ok

# UDIFF+ PDIFF + DIVUDIFF
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_u_diff...)
push!(wrc[2], r_u_diff...)
push!(wrc[3], c_u_diff...)
push!(wrc[1], w_p_diff...)
push!(wrc[2], r_p_diff...)
push!(wrc[3], c_p_diff...)
push!(wrc[1], w_divu_diff...)
push!(wrc[2], r_divu_diff...)
push!(wrc[3], c_divu_diff...)
A = assemble_matrix(assem, wrc)
res_stab_updivu_diff = compute_quantities(problem,A,b,dΩ) #ok

# UFULL + DIVUFULL + "MIX"
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_u_full...)
push!(wrc[2], r_u_full...)
push!(wrc[3], c_u_full...)
push!(wrc[1], w_divu_full...)
push!(wrc[2], r_divu_full...)
push!(wrc[3], c_divu_full...)
push!(wrc[1], w_divuq_full...)
push!(wrc[2], r_divuq_full...)
push!(wrc[3], c_divuq_full...)
push!(wrc[1], w_pdivv_full...)
push!(wrc[2], r_pdivv_full...)
push!(wrc[3], c_pdivv_full...)
A = assemble_matrix(assem, wrc)
res_stab_udivumix_full = compute_quantities(problem,A,b,dΩ)

# UDIFF + DIVUDIFF + "MIX"
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_u_diff...)
push!(wrc[2], r_u_diff...)
push!(wrc[3], c_u_diff...)
push!(wrc[1], w_divu_diff...)
push!(wrc[2], r_divu_diff...)
push!(wrc[3], c_divu_diff...)
A = assemble_matrix(assem, wrc)
res_stab_udivumix_rhs_diff= compute_quantities(problem,A,b_rhsq_diff,dΩ)
res_stab_udivumix_diff = compute_quantities(problem,A,b,dΩ) # this is the same for rhs_func = constant
@assert res_stab_udivumix_rhs_diff == res_stab_udivumix_diff

# UFULL + DIVUFULL + "MIX"
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))
push!(wrc[1], w_u_full...)
push!(wrc[2], r_u_full...)
push!(wrc[3], c_u_full...)
push!(wrc[1], w_divu_full...)
push!(wrc[2], r_divu_full...)
push!(wrc[3], c_divu_full...)
A = assemble_matrix(assem, wrc)
res_stab_udivumix_rhs_full= compute_quantities(problem,A,b_rhsq_full,dΩ)
res_stab_udivumix_full = compute_quantities(problem,A,b,dΩ) # this is the same for rhs_func = constant
@assert res_stab_udivumix_rhs_full == res_stab_udivumix_full

## GENERATE EXCEL SHEET RESULTS``
print("results")
res_nostab
res_stab_u_full
res_stab_u_diff
res_stab_p_full
res_stab_p_diff
res_stab_divu_full
res_stab_divu_diff
res_stab_pdivv_full # THIS IS OK. 
res_stab_divuq_full_rhs # ERRORS SEEM TO BE OF SAME ORDER NOW.

res_stab_updivu_full
res_stab_updivu_diff
res_stab_udivumix_full
res_stab_udivumix_rhs_diff
res_stab_udivumix_rhs_full