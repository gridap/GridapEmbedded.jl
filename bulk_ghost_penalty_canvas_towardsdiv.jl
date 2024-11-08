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

# Set up global spaces
Ωhact = Triangulation(cutdisk,ACTIVE)
degree=2*(order+1)
dΩhact = Measure(Ωhact,degree)

V = FESpace(Ωhact, ReferenceFE(raviart_thomas,Float64,order),conformity=:HDiv)
Q = FESpace(Ωhact, ReferenceFE(lagrangian,Float64,order), conformity=:L2)
U = TrialFESpace(V)
P = TrialFESpace(Q)
Y = MultiFieldFESpace([V, Q])
X = MultiFieldFESpace([U, P])
dx    = get_trial_fe_basis(X)
dy    = get_fe_basis(Y)
du,dp = dx
dv,dq = dy

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

## TESTING THE PROJECTION of a FE function


# writevtk(Ωhact,"output" , cellfields=["div_u"=>∇⋅uh,"pproj"=>testFEfunction])

function divergence_in_pressure_space_fe_function(uh::Gridap.FESpaces.SingleFieldFEFunction,P,dq)
    dq_data = Gridap.CellData.get_data(_get_single_field_fe_basis(dq))
    div_dv = ∇⋅(uh) ## We should be fine here.
    pdofs=Gridap.FESpaces.get_fe_dof_basis(P)

    # DoFs of the FEFunction cell-wise
    div_dv_pdofs_values_cell_wise=pdofs(div_dv)

    # Gather cell-wise DoFs into global DoFs vector
    div_dv_pdofs_values=Vector{Float64}(undef, num_free_dofs(P))
    Gridap.FESpaces.gather_free_values!(div_dv_pdofs_values,P,div_dv_pdofs_values_cell_wise)

    # free_dof_values is actually a global DoFs vector (!!!)
    FEFunction(P,div_dv_pdofs_values)
end
testxh=FEFunction(X, rand(num_free_dofs(X)))
testuh,testph = testxh
testFEfunction = divergence_in_pressure_space_fe_function(testuh,P, dq)
eh = testFEfunction-(∇⋅testuh)
norm_testFEfunction = sum(∫(eh⋅eh)dΩhact)

# TODO: the following functions should be tested and can then replace the lines below these commented lines.
# function divergence_in_pressure_space_basis(uh::FEBasis,P,dq)
#     dq_data = Gridap.CellData.get_data(_get_single_field_fe_basis(dq))
#     div_dv = ∇⋅(_get_single_field_fe_basis(uh))
#     pdofs=Gridap.FESpaces.get_fe_dof_basis(P)
#     div_dv_pdofs_values=pdofs(div_dv)
#     div_dv_in_pressure_space_cell_array=lazy_map(Gridap.Fields.linear_combination,
#                                                  div_dv_pdofs_values,
#                                                  dq_data)
#     if (Gridap.FESpaces.BasisStyle(uh)==Gridap.FESpaces.TrialBasis())
#        div_dv_in_pressure_space_cell_array=lazy_map(transpose,div_dv_in_pressure_space_cell_array)
#     end
#     div_dv_in_pressure_space_single_field=
#                                Gridap.FESpaces.SingleFieldFEBasis(div_dv_in_pressure_space_cell_array,
#                                                                  get_triangulation(P),
#                                                                  BasisStyle(uh),
#                                                                  Gridap.CellData.ReferenceDomain())
#     div_dv_in_pressure_space=
#       Gridap.MultiField.MultiFieldFEBasisComponent(div_dv_in_pressure_space_single_field,1,2)
# end
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
### Compute Π_Q_bb(div_du)
dqbb_Ωagg_cells =change_domain_bb_to_agg_cells(qbb,ref_agg_cell_to_ref_bb_map,Ωagg_cells,agg_cells_to_aggregate)

# OLD ROUTINE:
# agg_cells_rhs_contribs=get_array(∫(dqbb_Ωagg_cells⋅div_du_in_pressure_space_single_field)dΩagg_cells)
# ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(P_agg_cells_local_dof_ids,agg_cells_rhs_contribs) ## TODO: prep to assemble vectors
# rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)
test_agg_cells_rhs_contribs=get_array(∫(dqbb_Ωagg_cells⋅testFEfunction)dΩagg_cells)
test_ass_rhs_map=BulkGhostPenaltyAssembleRhsMap(P_agg_cells_local_dof_ids,test_agg_cells_rhs_contribs)
test_rhs=lazy_map(test_ass_rhs_map,aggregate_to_local_cells)

test_div_dv_l2_proj_bb_dofs=lazy_map(\,p_lhs,test_rhs)

# TODO: can be removed if we fix the BulkGhostPenaltyAssembleRhsMap.
# transposed_test_div_dv_l2_proj_bb_dofs= lazy_map(x->collect(reshape(x,(size(x,1)))),lazy_map(transpose,test_div_dv_l2_proj_bb_dofs))
# test_div_dv_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination,
#                                 transposed_test_div_dv_l2_proj_bb_dofs,
#                                 Gridap.CellData.get_data(qbb))

test_div_dv_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination,
                                test_div_dv_l2_proj_bb_dofs,
                                Gridap.CellData.get_data(qbb))
test_div_dv_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),
                                lazy_map(Reindex(test_div_dv_l2_proj_bb_array),agg_cells_to_aggregate),
                                ref_agg_cell_to_ref_bb_map)
div_uh_proj_bb=Gridap.CellData.GenericCellField(test_div_dv_l2_proj_bb_array_agg_cells,dΩagg_cells.quad.trian,ReferenceDomain())

x=get_cell_points(dΩagg_cells.quad.trian)
div_uh_proj_bb(x)[1]

res = testFEfunction - div_uh_proj_bb
norm = sum(∫((res)*(res))*dΩagg_cells)