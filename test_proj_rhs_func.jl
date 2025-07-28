using Gridap
using GridapEmbedded

include("./src/BGP/BGP.jl")

# Problem selection
problem = 0 # 0 = Manufactured solution (L2-like projection), 1 = Darcy problem 

# Manufactured solution
order  = 1
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

# # LHS of L2 projection on bounding boxes.

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

##==============================================================================##
#                                OBJECTIVE
##==============================================================================##
#=
  Compute and assemble the bulk penalty stabilization term

  γ ∫( (v)⊙(rhs_function - proj_rhs_function) )*dD

  with `v` the test basis functions (which may be components of a MultiField). The L2 projections `proj_rhs_function` can be computed via `setup_L2_proj_in_bb_space`. `dD` is the measure of the cut cells or full aggregate, that is `dΩbg_cut_cells` or `dΩbg_agg_cells`. The function `rhs_function` is the forcing term appearing on the right hand-side of the PDE. 

  Here, we are interested in using the pressure basis functions for the test space, thus v = dq. In addition, dD = dΩbg_cut_cells.
=#

# Manually set up the arrays that collect_cell_vector would return automatically
test_w = []
test_r = []

# Pick a rhs function
# rhs_func(x) = 10.0*x[1]^2 + sin(x[2])
rhs_func(x) = 10.0*x[1]
# rhs_func(x) = 1000.0

##==============================================================================##
#                   First term: γ ∫( (v)⊙(rhs_function) )*dD
##==============================================================================##
test_rhs_vec_contribs1=get_array(∫(γ*(rhs_func⊙dq))*dΩbg_cut_cells) # 16-element array, one entry per T ∈ Ωbg_cut_cells, with 1-element Vector (for pressure) per cell T.
push!(test_w, test_rhs_vec_contribs1)

if (_is_multifield_fe_basis_component(dq))
   nfields=_nfields(dq)
   fieldid=_fieldid(dq)
   testP_Ωbg_cut_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),P_Ωbg_cut_cell_dof_ids)
end
push!(test_r, testP_Ωbg_cut_cell_dof_ids) 

##==============================================================================##
#                   Second term: - γ ∫( (v)⊙(proj_rhs_function) )*dD
##==============================================================================##

## (1) Set-up proj_rhs_function using 
#=

returns the L2 projection in the bounding box space of the rhs function. That is, the L2 projection is defined through

        ∫ (rhs_func - Π_{Zbb}(operation(u)) zbb dΩbg_agg_cells = 0   ∀ zbb ∈ Zbb(T), 
        
with T ∈ T_agg and `Zbb` the bounding box space. Note that Ωbg_agg_cells ⊆ Ωbg_bb. Additionally, Π_{Zbb}(operation(u)) appears in the lhs as the trial function wbb ∈ Zbb.

=#
test_wbb     = pbb      # Bounding box trial space (for pressure)
test_zbb     = qbb      # Bounding box test space (for pressure)   

test_wbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(test_wbb, ref_agg_cell_to_ref_bb_map, Ωbg_agg_cells, agg_cells_to_aggregate)   # trial
test_zbb_Ωbg_agg_cells=change_domain_bb_to_agg_cells(test_zbb, ref_agg_cell_to_ref_bb_map, Ωbg_agg_cells, agg_cells_to_aggregate)   # test
test_agg_cells_to_lhs_contribs=get_array(∫(test_zbb_Ωbg_agg_cells⋅test_wbb_Ωbg_agg_cells)dΩbg_agg_cells)
test_ass_lhs_map= BulkGhostPenaltyAssembleLhsMap(test_agg_cells_to_lhs_contribs)
test_lhs        = lazy_map(test_ass_lhs_map,aggregate_to_local_cells)
test_agg_cells_rhs_contribs=get_array(∫(test_zbb_Ωbg_agg_cells⋅rhs_func)dΩbg_agg_cells)
test_rhs = lazy_map(sum, 
                    lazy_map(Broadcasting(Reindex(test_agg_cells_rhs_contribs)), 
                    aggregate_to_local_cells))          
test_f_proj_Zbb_dofs=lazy_map(\,test_lhs,test_rhs)
test_f_proj_Zbb_array=lazy_map(Gridap.Fields.linear_combination, 
                                  test_f_proj_Zbb_dofs, 
                               Gridap.CellData.get_data(test_zbb))

# Change domain of proj_op_u and proj_op_v from Ωbb to Ωbg_agg_cells
test_f_proj_Zbb_array_agg_cells=lazy_map(Broadcasting(∘), 
                          lazy_map(Reindex(test_f_proj_Zbb_array),agg_cells_to_aggregate), ref_agg_cell_to_ref_bb_map)
test_f_proj_Zbb_agg_cells = 
    Gridap.CellData.GenericCellField(test_f_proj_Zbb_array_agg_cells,Ωbg_agg_cells,ReferenceDomain())

## (2) Define γ ∫( (v)⊙(proj_rhs_function) )*dD
dq_on_Ωbg_agg_cells = Gridap.CellData.change_domain(dq,Ωbg_agg_cells,ReferenceDomain())
test_rhs_vec_contribs2 = get_array(∫((-1.0)*γ*(test_f_proj_Zbb_agg_cells)⋅dq_on_Ωbg_agg_cells)*dΩbg_cut_cells)
push!(test_w, test_rhs_vec_contribs2)
push!(test_r, testP_Ωbg_cut_cell_dof_ids) # same dofs as as contr1? TOCHECK!

# ### TEST FOR INDIVIDUAL CONTRIBUTIONS: I expect these to cancel out
# for i= 1:16
#    println("===== i: $i")
#    println("contr1: $(test_rhs_vec_contribs1[i][2])")
#    println("contr2: $(test_rhs_vec_contribs2[i][2])")
# end

# ### TEST FOR SUM OF CONTRIBUTIONS: I expect these to cancel out
# ([sum(get_array(∫(test_f_proj_Zbb_agg_cells)dΩbg_cut_cells)[i]) for i=1:16])
# sum_cut_cells_rhs = sum([sum(get_array(∫(test_f_proj_Zbb_agg_cells)dΩbg_cut_cells)[i]) for i=1:16])
# sum((get_array(∫(rhs_func)*dΩbg_cut_cells)))
# difference = abs(sum_cut_cells_rhs - sum((get_array(∫(rhs_func)*dΩbg_cut_cells)))) 

#### -- TEST PROBLEM -- ####
test_a((u,p),(v,q))=∫(u⋅v+q*p)dΩ 
uex(x) = -VectorValue(2.0,2.0)
test_l((v,q))=∫(uex⋅v+0.0⋅q)dΩ

test_mat_A=Gridap.FESpaces.collect_cell_matrix(X,Y,test_a(dx,dy))
test_vec_b=Gridap.FESpaces.collect_cell_vector(Y,test_l(dy))
push!(test_vec_b[1], test_w...)
push!(test_vec_b[2], test_r...)

assem=SparseMatrixAssembler(X,Y)

test_A = assemble_matrix(assem, test_mat_A)
test_b = assemble_vector(assem, test_vec_b)

test_sol = test_A\test_b

test_sol_FE = FEFunction(X, test_sol)
test_sol_u, test_sol_p = test_sol_FE

norm_sol_u = (sum(∫(test_sol_u⋅test_sol_u)*dΩ))
norm_sol_p = (sum(∫(test_sol_p*test_sol_p)*dΩ))

err_sol_u = uex - test_sol_u
err_sol_p = 0.0 - test_sol_p
L2error_sol_u = (sum(∫(err_sol_u⋅err_sol_u)*dΩ))
L2error_sol_p = (sum(∫(err_sol_p⋅err_sol_p)*dΩ))

##################### TRY TO DEFINE FUNCTIONS TO COLLECT VEC_STAB
proj_rhs_func = setup_L2_proj_in_bb_space(dΩbg_agg_cells,
   ref_agg_cell_to_ref_bb_map, # map
   agg_cells_to_aggregate,     # 
   aggregate_to_local_cells,   # 
   rhs_func,                   # Function to be projected 
   pbb,                        # Trial basis of bounding box space Zbb
   qbb) 
# #THIS TEST SEEMS TO IMPLY THAT ALL IS FINE.
# for i=1:16
#       println(get_array(∫(proj_rhs_func)dΩbg_cut_cells)[i])
#       println(get_array(∫(proj_rhs_func)dΩbg_cut_cells)[i]-get_array(∫(test_f_proj_Zbb_agg_cells)dΩbg_cut_cells)[i])
# end

wvec, rvec =  bulk_ghost_penalty_stabilization_collect_cell_vector_on_D(dΩbg_cut_cells,
                                          γ,  
                                          dq, 
                                          testP_Ωbg_cut_cell_dof_ids,
                                          identity,
                                          rhs_func,
                                          proj_rhs_func)

#### -- TEST PROBLEM -- #### MET wvec
test_a((u,p),(v,q))=∫(u⋅v+q*p)dΩ 
uex(x) = -VectorValue(2.0,2.0)
test_l((v,q))=∫(uex⋅v+0.0⋅q)dΩ

test_mat_A=Gridap.FESpaces.collect_cell_matrix(X,Y,test_a(dx,dy))
test_vec_b=Gridap.FESpaces.collect_cell_vector(Y,test_l(dy))
push!(test_vec_b[1], wvec...)
push!(test_vec_b[2], rvec...)

assem=SparseMatrixAssembler(X,Y)

test_A = assemble_matrix(assem, test_mat_A)
test_b = assemble_vector(assem, test_vec_b)

test_sol = test_A\test_b

test_sol_FE = FEFunction(X, test_sol)
test_sol_u, test_sol_p = test_sol_FE

norm_sol_u = (sum(∫(test_sol_u⋅test_sol_u)*dΩ))
norm_sol_p = (sum(∫(test_sol_p*test_sol_p)*dΩ))

err_sol_u = uex - test_sol_u
err_sol_p = 0.0 - test_sol_p
L2error_sol_u = (sum(∫(err_sol_u⋅err_sol_u)*dΩ))
L2error_sol_p = (sum(∫(err_sol_p⋅err_sol_p)*dΩ))

####
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

if (_is_multifield_fe_basis_component(dp))
   nfields=_nfields(dp)
   fieldid=_fieldid(dp)
   P_cut_cells_to_aggregate_dof_ids=
   lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),P_cut_cells_to_aggregate_dof_ids)
end
wvec2, rvec2 =  bulk_ghost_penalty_stabilization_collect_cell_vector_on_D(dΩbg_cut_cells,
                                          γ,  
                                          dq,
                                          dq_proj_Qbb,
                                          testP_Ωbg_cut_cell_dof_ids,
                                          P_cut_cells_to_aggregate_dof_ids,
                                          identity,
                                          rhs_func,
                                          proj_rhs_func)

#### -- TEST PROBLEM -- #### MET wvec2
test_a((u,p),(v,q))=∫(u⋅v+q*p)dΩ 
uex(x) = -VectorValue(2.0,2.0)
test_l((v,q))=∫(uex⋅v+0.0⋅q)dΩ

test_mat_A=Gridap.FESpaces.collect_cell_matrix(X,Y,test_a(dx,dy))
test_vec_b=Gridap.FESpaces.collect_cell_vector(Y,test_l(dy))
push!(test_vec_b[1], wvec2...)
push!(test_vec_b[2], rvec2...)

assem=SparseMatrixAssembler(X,Y)

test_A = assemble_matrix(assem, test_mat_A)
test_b = assemble_vector(assem, test_vec_b)

test_sol = test_A\test_b

test_sol_FE = FEFunction(X, test_sol)
test_sol_u, test_sol_p = test_sol_FE

norm_sol_u = (sum(∫(test_sol_u⋅test_sol_u)*dΩ))
norm_sol_p = (sum(∫(test_sol_p*test_sol_p)*dΩ))

err_sol_u = uex - test_sol_u
err_sol_p = 0.0 - test_sol_p
L2error_sol_u = (sum(∫(err_sol_u⋅err_sol_u)*dΩ))
L2error_sol_p = (sum(∫(err_sol_p⋅err_sol_p)*dΩ))