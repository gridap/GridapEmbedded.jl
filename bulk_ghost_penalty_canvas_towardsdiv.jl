using Gridap
using GridapEmbedded
using FillArrays
using LinearAlgebra

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

"""
  Creates an array of arrays with as many entries 
  as aggregates. For each aggregate, the array 
  contains the global cell IDs of that cells in the background 
  model that belong to the same aggregate

  TO-DO: with efficiency in mind we may want to store this 
         array of arrays as a Gridap.Arrays.Table.
"""
function setup_aggregate_to_cells(aggregates)
  
  size_aggregates=Dict{Int,Int}()
  for (i,agg) in enumerate(aggregates)
    if agg>0
        if !haskey(size_aggregates,agg)
            size_aggregates[agg]=1
        else 
            size_aggregates[agg]+=1
        end
    end
  end   

  touched=Dict{Int,Int}()
  aggregate_to_cells=Vector{Vector{Int}}()
  current_aggregate=1
  for (i,agg) in enumerate(aggregates)
    if agg>0
        if (size_aggregates[agg]>1)
            if !haskey(touched,agg)
                push!(aggregate_to_cells,[i])
                touched[agg]=current_aggregate
                current_aggregate+=1
            else 
                push!(aggregate_to_cells[touched[agg]],i)
            end
        end 
    end
  end
  aggregate_to_cells
end 

function setup_aggregates_bounding_box_model(bgmodel, aggregate_to_cells)
    g=get_grid(bgmodel)
    cell_coords=get_cell_coordinates(g)
    D=num_dims(bgmodel)
    xmin=Vector{Float64}(undef,D)
    xmax=Vector{Float64}(undef,D)
    
    # Compute coordinates of the nodes defining the bounding boxes
    bounding_box_node_coords=
    Vector{Point{D,Float64}}(undef,length(aggregate_to_cells)*2^D)
    ptr  = [ (((i-1)*2^D)+1) for i in 1:length(aggregate_to_cells)+1 ]
    data = collect(1:length(bounding_box_node_coords))
    bounding_box_node_ids = Gridap.Arrays.Table(data,ptr)
    for (agg,cells) in enumerate(aggregate_to_cells)
        p=first(cell_coords[cells[1]])
        for i in 1:D
            xmin[i]=p[i]
            xmax[i]=p[i]
        end
        for cell in cells
            for p in cell_coords[cell]
                for i in 1:D
                    xmin[i]=min(xmin[i],p[i])
                    xmax[i]=max(xmax[i],p[i])
                end 
            end 
        end
        bounds         = [(xmin[i], xmax[i]) for i in 1:D]
        point_iterator = Iterators.product(bounds...)
        bounding_box_node_coords[bounding_box_node_ids[agg]] = 
                reshape([Point(p...) for p in point_iterator],2^D)
    end

    # Set up the discrete model of bounding boxes
    HEX_AXIS=1
    polytope=Polytope(Fill(HEX_AXIS,D)...)
    scalar_reffe=ReferenceFE(polytope,lagrangian,Float64,1)
    cell_types=fill(1,length(bounding_box_node_ids))
    cell_reffes=[scalar_reffe]
    grid = Gridap.Geometry.UnstructuredGrid(bounding_box_node_coords,
                                            bounding_box_node_ids,
                                            cell_reffes,
                                            cell_types,
                                            Gridap.Geometry.Oriented())
    Gridap.Geometry.UnstructuredDiscreteModel(grid)
end 

aggregate_to_cells=setup_aggregate_to_cells(aggregates)
aggregates_bounding_box_model=
       setup_aggregates_bounding_box_model(bgmodel,aggregate_to_cells)

# colors = color_aggregates(aggregates,bgmodel)
# writevtk(Triangulation(bgmodel),"trian",celldata=["cellin"=>aggregates,"color"=>colors])       
# writevtk(aggregates_bounding_box_model, "bb_model")
# writevtk(bgmodel, "bg_model")

"""
  Changes the domain of a trial/test basis defined on 
  the reference space of bounding boxes to the reference 
  space of the agg cells

  TO-DO: in the future, for system of PDEs (MultiField) we should 
         also take care of blocks (BlockMap)
"""
function change_domain_bb_to_agg_cells(basis_bb, 
                                       ref_agg_cell_to_ref_bb_map, 
                                       Ωagg_cells,
                                       agg_cells_to_aggregate)
    @assert num_cells(Ωagg_cells)==length(ref_agg_cell_to_ref_bb_map)
    @assert Gridap.CellData.DomainStyle(basis_bb)==ReferenceDomain()
    bb_basis_style = Gridap.FESpaces.BasisStyle(basis_bb)  
    bb_basis_array = Gridap.CellData.get_data(basis_bb)
    if (bb_basis_style==Gridap.FESpaces.TrialBasis())
        # Remove transpose map; we will add it later
        @assert isa(bb_basis_array,Gridap.Arrays.LazyArray)
        @assert isa(bb_basis_array.maps,Fill)
        @assert isa(bb_basis_array.maps.value,typeof(transpose))
        bb_basis_array=bb_basis_array.args[1]
    end
    
    bb_basis_array_to_Ωagg_cells_array = lazy_map(Reindex(bb_basis_array),agg_cells_to_aggregate)
    bb_basis_array_to_Ωagg_cells_array = lazy_map(Broadcasting(∘),
                                                  bb_basis_array_to_Ωagg_cells_array,
                                                  ref_agg_cell_to_ref_bb_map)
    if (bb_basis_style==Gridap.FESpaces.TrialBasis())
        # Add transpose 
        bb_basis_array_to_Ωagg_cells_array=lazy_map(transpose, bb_basis_array_to_Ωagg_cells_array)
    end    

    Gridap.CellData.GenericCellField(bb_basis_array_to_Ωagg_cells_array,
                                     Ωagg_cells,
                                     ReferenceDomain())
end


# Set up objects required to compute both LHS and RHS of the L2 projection

# Set up global spaces 
Ωhact = Triangulation(cutdisk,ACTIVE)

V = FESpace(Ωhact, ReferenceFE(raviart_thomas,Float64,order),conformity=:HDiv)
Q = FESpace(Ωhact, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
U = TrialFESpace(V)
P = TrialFESpace(Q)
Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])

# Generate an array with the global IDs of the cells 
# that belong to an aggregrate. From now on, we will 
# use the terminology "agg_cells" to refer to those 
# cells of the background model that belong to an aggregate 
# (i.e., they can be either cut or interior cells)
agg_cells=Vector{Int}()
for cells in aggregate_to_cells
    append!(agg_cells,cells)
end

# ref_agg_cell_to_agg_cell_map: \hat{K} -> K 
Ω=Triangulation(bgmodel)
Ωagg_cells=view(Ω,agg_cells)
ref_agg_cell_to_agg_cell_map=get_cell_map(Ωagg_cells)

# Generate an array that given the local ID of an "agg_cell"
# returns the ID of the aggregate to which it belongs
# (i.e., flattened version of aggregate_to_cells)
agg_cells_to_aggregate=Vector{Int}()
for (i,cells) in enumerate(aggregate_to_cells)
    for _ in cells 
        push!(agg_cells_to_aggregate,i)
    end
end 

# ref_agg_cell_to_ref_bb_map: \hat{K} -> K -> bb -> \hat{bb}
bb_to_ref_bb=lazy_map(Gridap.Fields.inverse_map,get_cell_map(aggregates_bounding_box_model))
bb_to_ref_bb_agg_cells=lazy_map(Reindex(bb_to_ref_bb),agg_cells_to_aggregate)  
ref_agg_cell_to_ref_bb_map=
    lazy_map(Broadcasting(∘),bb_to_ref_bb_agg_cells,ref_agg_cell_to_agg_cell_map)

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


aggregate_to_local_cells=copy(aggregate_to_cells)
current_local_cell=1
for (i,cells) in enumerate(aggregate_to_local_cells)
    for j in 1:length(cells)
        cells[j]=current_local_cell
        current_local_cell+=1
    end
end

function set_up_bulk_ghost_penalty_lhs(aggregates_bounding_box_model,
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

function _restrict_to_block(cell_dof_ids::Gridap.Arrays.LazyArray{<:Fill{<:Gridap.Fields.BlockMap}}, blockid)
    map=cell_dof_ids.maps.value 
    @assert length(map.size)==1
    @assert blockid >= 1
    @assert blockid <= map.size[1]
    cell_dof_ids.args[blockid]
end 

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

function compute_agg_cells_local_dof_ids(agg_cells_dof_ids, aggregate_to_agg_cells)
    agg_cells_local_dof_ids=copy(agg_cells_dof_ids)
    current_cell=1
    for agg_cells in aggregate_to_agg_cells
        g2l=Dict{Int32,Int32}()
        current_local_dof=1
        for (i,_) in enumerate(agg_cells)
            current_cell_dof_ids=agg_cells_dof_ids[current_cell]
            for (j, dof) in enumerate(current_cell_dof_ids)
                if !(dof in keys(g2l))
                    g2l[dof]=current_local_dof
                    agg_cells_local_dof_ids[current_cell][j]=current_local_dof
                    current_local_dof+=1
                else 
                    agg_cells_local_dof_ids[current_cell][j]=g2l[dof]
                end 
            end 
            current_cell+=1
            println(agg_cells_local_dof_ids)
        end
    end
    agg_cells_local_dof_ids
end

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


function compute_aggregate_dof_ids(agg_cells_dof_ids, aggregate_to_agg_cells)
    aggregate_dof_ids=Vector{Vector{Int}}(undef, length(aggregate_to_agg_cells))
    current_aggregate=1
    current_cell=1
    for agg_cells in aggregate_to_agg_cells
        current_aggregate_dof_ids=Int[]
        for (i,_) in enumerate(agg_cells)
            current_cell_dof_ids=agg_cells_dof_ids[current_cell]
            for (j, dof) in enumerate(current_cell_dof_ids)
                if !(dof in current_aggregate_dof_ids)
                    push!(current_aggregate_dof_ids, dof)
                end 
            end 
            current_cell+=1
        end 
        aggregate_dof_ids[current_aggregate]=current_aggregate_dof_ids
        current_aggregate+=1
    end
    aggregate_dof_ids
end

function _get_single_field_fe_basis(a::Gridap.MultiField.MultiFieldFEBasisComponent)
    a.single_field
end
function _get_single_field_fe_basis(a)
    a
end
function _is_multifield_fe_basis_component(a::Gridap.MultiField.MultiFieldFEBasisComponent)
    true
end
function _is_multifield_fe_basis_component(a)
    false
end
function _nfields(a::Gridap.MultiField.MultiFieldFEBasisComponent)
    a.nfields
end
function _fieldid(a::Gridap.MultiField.MultiFieldFEBasisComponent)
    a.fieldid
end

"""
  dv, du: Test and trial basis functions. They may be components of a MultiFieldCellField

  # Compute and assemble the bulk penalty stabilization term
  # ∫( (dv-dv_l2_proj_agg_cells)*(du-du_l2_proj_agg_cells))*dΩ_agg_cells (long version)
  # ∫( (dv)*(du-du_l2_proj_agg_cells))*dΩ_agg_cells (simplified, equivalent version)
"""
function interior_bulk_penalty_stabilization_collect_cell_matrix(agg_cells_to_aggregate,
                                                                 ref_agg_cell_to_ref_bb_map,
                                                                 dΩagg_cells,
                                                                 dv,    # Test basis
                                                                 du,    # Trial basis (to project)
                                                                 dvbb,  # Bounding box space test basis
                                                                 lhs,
                                                                 Ωagg_cell_dof_ids,
                                                                 agg_cells_local_dof_ids,
                                                                 agg_cells_to_aggregate_dof_ids,
                                                                 h_U,
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
  
U_aggregate_dof_ids=compute_aggregate_dof_ids(U_Ωagg_cell_dof_ids,aggregate_to_cells)
U_agg_cells_to_aggregate_dof_ids=lazy_map(Reindex(U_aggregate_dof_ids),agg_cells_to_aggregate)

P_aggregate_dof_ids=compute_aggregate_dof_ids(P_Ωagg_cell_dof_ids,aggregate_to_cells)
P_agg_cells_to_aggregate_dof_ids=lazy_map(Reindex(P_aggregate_dof_ids),agg_cells_to_aggregate)

γ = 10.0 # Interior bulk-penalty stabilization parameter
         # (@amartinhuertas no idea what a reasonable value is)

h_U = set_up_h_U(aggregates_bounding_box_model, agg_cells_to_aggregate, Ωagg_cells)         

# Manually set up the arrays that collect_cell_matrix would return automatically
wp,rp,cp=interior_bulk_penalty_stabilization_collect_cell_matrix(agg_cells_to_aggregate,
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

wu,ru,cu=interior_bulk_penalty_stabilization_collect_cell_matrix(agg_cells_to_aggregate,
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

Ωcut = Triangulation(cutdisk,PHYSICAL)
dΩcut = Measure(Ωcut,degree)
# # Set up global projection matrix 
# Ωcut = Triangulation(cutdisk,PHYSICAL)
# dΩcut = Measure(Ωcut,degree)
# a_nodiv((u,p),(v,q))=∫(v⋅u+q*p)dΩcut 
# wrc_nodiv=Gridap.FESpaces.collect_cell_matrix(X,Y,a_nodiv(dx,dy))

# assem_nodiv=SparseMatrixAssembler(X,Y)
# Anostab_nodiv=assemble_matrix(assem_nodiv, wrc_nodiv)
# cond(Array(Anostab_nodiv))

# # Add the bulk penalty stabilization term to wrc_nodiv
# push!(wrc_nodiv[1], wp...)
# push!(wrc_nodiv[2], rp...)
# push!(wrc_nodiv[3], cp...)

# push!(wrc_nodiv[1], wu...)
# push!(wrc_nodiv[2], ru...)
# push!(wrc_nodiv[3], cu...)

# Awithstab_nodiv=assemble_matrix(assem_nodiv, wrc_nodiv)
# cond(Array(Awithstab_nodiv))

# # Set up rhs global projection 
# l_nodiv((v,q))=∫(v⋅uex+q*pex)dΩcut
# b_nodiv = assemble_vector(l_nodiv, Y)

# global_l2_proj_dofs_nodiv = Awithstab_nodiv\b_nodiv
# xh_nodiv = FEFunction(X, global_l2_proj_dofs_nodiv)
# uh_nodiv,ph_nodiv = xh_nodiv

# euh_nodiv = uex-uh_nodiv
# # @assert sum(∫(euh_nodiv⋅euh_nodiv)*dΩcut) < 1.0e-12

# eph_nodiv = pex-ph_nodiv
# # @assert sum(∫(eph_nodiv*eph_nodiv)*dΩcut) < 1.0e-12
# eph_nodiv_norm   = sum(∫(eph_nodiv*eph_nodiv)*dΩcut)
# euh_nodiv_norm   = sum(∫(euh_nodiv⋅euh_nodiv)*dΩcut)
# # k=0 yields eph = 0.00012 and euh = 5.5e-30
# edivuh_nodiv = divuex-(∇⋅uh_nodiv)
# edivuh_nodiv_norm   = sum(∫(edivuh_nodiv⋅edivuh_nodiv)*dΩcut)

## TEST 1
a((u,p),(v,q))=∫(v⋅u+q*p+(∇⋅u)*(∇⋅v))dΩcut 
wrc=Gridap.FESpaces.collect_cell_matrix(X,Y,a(dx,dy))

assem=SparseMatrixAssembler(X,Y)
Anostab=assemble_matrix(assem, wrc)
cond_nostab = cond(Array(Anostab))

# Set up rhs global projection 
l((v,q))=∫(v⋅uex+q*pex+(∇⋅v)*divuex)dΩcut
b = assemble_vector(l, Y)

## NO STABILIZATION 
if cond_nostab<1e35
    global_l2_proj_dofs_nostab = Anostab\b
    xh_nostab = FEFunction(X, global_l2_proj_dofs_nostab)
    uh_nostab,ph_nostab = xh_nostab
    euh_nostab = uex-uh_nostab
    eph_nostab = pex-ph_nostab
    edivuh_nostab = divuex-(∇⋅uh_nostab)
    norm_euh_nostab = sum(∫(euh_nostab⋅euh_nostab)*dΩcut)
    norm_eph_nostab = sum(∫(eph_nostab*eph_nostab)*dΩcut)
    edivuh_nostab   = sum(∫(edivuh_nostab⋅edivuh_nostab)*dΩcut)
else
    println("Singular exception is raised as error when solving the system (as cond= $cond_nostab is high?)")
end
# k=0 yields euh = 1.2e-32, eph = 0.00017, edivuh = 6.8e-30, κ= 3.15e32
# k=0 yields euh = ?, eph = ?, edivuh = ?, κ= 4.3e36

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
# k=0 yields euh = 6.3e-31, eph = 0.00055, edivuh = 2.4e-28, κ= 11492
# k=2 yields euh = 5.3e-22, eph = 3.0e-18, edivuh = 3.0e-18, κ= 1.7e11

## +DIVU STABILIZATION 
## Compute ∫( (div_dv)*(div_du-div_du_l2_proj_agg_cells))*dΩ_agg_cells (simplified, equivalent version)
# Manually set up the arrays that collect_cell_matrix would return automatically
w_divu = []
r_divu = []
c_divu = []

## (I) Let's set up the div_dv_div_du_mat_contribs (first part of the integral, that is ∫( (div_dv)*(div_du)*dΩ_agg_cells
div_dv_div_du_mat_contribs=get_array(∫(γ*(∇⋅dv)⋅(∇⋅du))*dΩagg_cells)
U_Ωagg_cell_dof_ids = _restrict_to_block(Ωagg_cell_dof_ids, 1) 
if (_is_multifield_fe_basis_component(dv))
    nfields=_nfields(dv)
    fieldid=_fieldid(dv)
    U_Ωagg_cell_dof_ids=lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),U_Ωagg_cell_dof_ids)
end

push!(w_divu, div_dv_div_du_mat_contribs)
push!(r_divu, U_Ωagg_cell_dof_ids)
push!(c_divu, U_Ωagg_cell_dof_ids)

## Set up  ∫( (div_dv)*(-div_du_l2_proj_agg_cells))*dΩ_agg_cells
div_dv=∇⋅(_get_single_field_fe_basis(dv))
pdofs=Gridap.FESpaces.get_fe_dof_basis(P)
div_dv_pdofs_values=pdofs(div_dv)
div_dv_in_pressure_space_cell_array=lazy_map(Gridap.Fields.linear_combination,
                           div_dv_pdofs_values,
                           Gridap.CellData.get_data(_get_single_field_fe_basis(dq)))
div_dv_in_pressure_space_single_field= Gridap.FESpaces.SingleFieldFEBasis(div_dv_in_pressure_space_cell_array,
                                              get_triangulation(P),
                                              Gridap.FESpaces.TestBasis(),
                                              Gridap.CellData.ReferenceDomain())
div_dv_in_pressure_space=Gridap.MultiField.MultiFieldFEBasisComponent(div_dv_in_pressure_space_single_field,1,2)                                          

div_du_in_pressure_space_cell_array=lazy_map(transpose,div_dv_in_pressure_space_cell_array)
div_du_in_pressure_space_single_field=Gridap.FESpaces.SingleFieldFEBasis(div_du_in_pressure_space_cell_array,
                                                                 get_triangulation(P),
                                                                 Gridap.FESpaces.TrialBasis(),
                                                                 Gridap.CellData.ReferenceDomain())
div_du_in_pressure_space=Gridap.MultiField.MultiFieldFEBasisComponent(div_du_in_pressure_space_single_field,1,2)

div_dv_div_du_in_pressure_space_mat_contribs=get_array(∫(γ*div_dv_in_pressure_space⋅div_du_in_pressure_space)*dΩagg_cells) #not sure if neccesary?
div_dv_div_du_mat_contribs ≈ div_dv_div_du_in_pressure_space_mat_contribs

### Compute Π_Q_bb(div_du)
dqbb_Ωagg_cells =change_domain_bb_to_agg_cells(qbb,ref_agg_cell_to_ref_bb_map,Ωagg_cells,agg_cells_to_aggregate)  
agg_cells_rhs_contribs=get_array(∫(dqbb_Ωagg_cells⋅div_du_in_pressure_space_single_field)dΩagg_cells) 
# agg_cells_rhs_contribs[1]
ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(P_agg_cells_local_dof_ids,agg_cells_rhs_contribs)
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

if (_is_multifield_fe_basis_component(dp))
    @assert _is_multifield_fe_basis_component(dq)
    @assert _nfields(dp)==_nfields(dq)
    nfields=_nfields(dp)
    fieldid=_fieldid(dp)
    div_du_l2_proj_bb_array_agg_cells=lazy_map(
                                Gridap.Fields.BlockMap((1,nfields),fieldid),
                                div_du_l2_proj_bb_array_agg_cells)
end 

div_du_l2_proj_agg_cells = Gridap.CellData.GenericCellField(div_du_l2_proj_bb_array_agg_cells, dΩagg_cells.quad.trian,ReferenceDomain())

# In the MultiField case, I have had to add this change domain 
# call before setting up the term right below. Otherwise, we get an error 
# when trying to multiply the fields. Not sure why this is happening
# dv_Ωagg_cells=Gridap.CellData.change_domain(dv,Ωagg_cells,ReferenceDomain())
div_dv_Ωagg_cells=Gridap.CellData.change_domain(∇⋅dv,Ωagg_cells,ReferenceDomain()) 

proj_div_dv_div_du_mat_contribs=get_array(∫(γ*(-1.0)*div_dv_Ωagg_cells*(div_du_l2_proj_agg_cells))*dΩagg_cells)

P_Ωagg_cell_dof_ids = _restrict_to_block(Ωagg_cell_dof_ids, 2)
P_aggregate_dof_ids=compute_aggregate_dof_ids(P_Ωagg_cell_dof_ids,aggregate_to_cells)
P_agg_cells_to_aggregate_dof_ids=lazy_map(Reindex(P_aggregate_dof_ids),agg_cells_to_aggregate)

if (_is_multifield_fe_basis_component(dp))
    nfields=_nfields(dp)
    fieldid=_fieldid(dp)
    P_agg_cells_to_aggregate_dof_ids=
        lazy_map(Gridap.Fields.BlockMap(nfields,fieldid),P_agg_cells_to_aggregate_dof_ids)
end

push!(w_divu, proj_div_dv_div_du_mat_contribs)
push!(r_divu, U_Ωagg_cell_dof_ids)              # test
push!(c_divu, P_agg_cells_to_aggregate_dof_ids) # trial 

## Now add div stabilization to A
push!(wrc[1], w_divu...)
push!(wrc[2], r_divu...)
push!(wrc[3], c_divu...)

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
# k=0 yields euh = 0.0014, eph = 0.47, edivuh = 0.47, κ= 8.51e11
# k=2 yields euh = 0.00013, eph = 0.088, edivuh = 0.088, κ= 2.9e14