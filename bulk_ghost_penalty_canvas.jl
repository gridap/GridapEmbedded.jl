using Gridap
using GridapEmbedded
using FillArrays
using LinearAlgebra

# Manufactured solution
order = 1
uex(x) = x[1]^order + x[2]^order

# Select geometry
R = 0.2
geom = disk(R, x0=Point(0.5,0.5))

# Setup background model
n=20
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
    current=1
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
        # TO-DO: this nested loop is only valid for D=2.
        #        Ideally, we should be able to generalize it for 
        #        arbitrary D
        for i in (xmin[2],xmax[2])
            for j in (xmin[1],xmax[1])
                bounding_box_node_coords[current]=Point(j,i)
                current+=1
            end 
        end
    end

    # Set up the discrete model of bounding boxes 
    ptr  = [ (((i-1)*2^D)+1) for i in 1:length(aggregate_to_cells)+1 ]
    data = collect(1:length(bounding_box_node_coords))
    bounding_box_node_ids = Gridap.Arrays.Table(data,ptr)

    # TO-DO: this polytope is only valid for D=2
    polytope=QUAD
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

colors = color_aggregates(aggregates,bgmodel)
writevtk(Triangulation(bgmodel),"trian",celldata=["cellin"=>aggregates,"color"=>colors])       
writevtk(aggregates_bounding_box_model, "bb_model")
writevtk(bgmodel, "bg_model")

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
reffe = ReferenceFE(lagrangian,Float64,order)
Vstd  = TestFESpace(Ωhact,reffe,conformity=:L2)
Ustd  = TrialFESpace(Vstd)

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
reffe=ReferenceFE(lagrangian,Float64,order) # Here we MUST use a Q space (not a P space!)
Vbb=FESpace(aggregates_bounding_box_model,reffe,conformity=:L2) # We need a DG space to represent the L2 projection
Ubb=TrialFESpace(Vbb)
ubb=get_trial_fe_basis(Ubb)
vbb=get_fe_basis(Vbb)

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
agg_cells_to_lhs_contribs=get_array(∫(vbb_Ωagg_cells*ubb_Ωagg_cells)dΩagg_cells)

aggregate_to_local_cells=copy(aggregate_to_cells)
current_local_cell=1
for (i,cells) in enumerate(aggregate_to_local_cells)
    for j in 1:length(cells)
        cells[j]=current_local_cell
        current_local_cell+=1
    end
end

# TO-DO: Better name?
struct AssembleLhsMap{A} <: Gridap.Fields.Map
    agg_cells_lhs_contribs::A
end

function _get_rank(::Type{Array{T,N}}) where {T,N}
  N
end

function Gridap.Fields.return_cache(m::AssembleLhsMap,cells)
    cache_unassembled_lhs=array_cache(m.agg_cells_lhs_contribs)
    T=eltype(m.agg_cells_lhs_contribs)
    evaluate_result=Gridap.Arrays.CachedArray(eltype(T),_get_rank(T))
    cache_unassembled_lhs,evaluate_result
end
  
function Gridap.Fields.evaluate!(cache,m::AssembleLhsMap,cells)
    cache_unassembled_lhs,result=cache
    contrib = getindex!(cache_unassembled_lhs,m.agg_cells_lhs_contribs,1)

    Gridap.Arrays.setsize!(result,size(contrib))
    result.array .= 0.0
    for (i,cell) in enumerate(cells)
        contrib = getindex!(cache_unassembled_lhs,m.agg_cells_lhs_contribs,cell)
        result.array .+= contrib
    end
    result.array
end

# Finally assemble LHS contributions
ass_lhs_map=AssembleLhsMap(agg_cells_to_lhs_contribs)
lhs=lazy_map(ass_lhs_map,aggregate_to_local_cells)

# Compute contributions to the RHS of the L2 projection
du    = get_trial_fe_basis(Ustd)
Ωagg_cell_dof_ids = get_cell_dof_ids(Ustd,Ωagg_cells)
agg_cells_rhs_contribs=get_array(∫(vbb_Ωagg_cells*du)dΩagg_cells)

# TO-DO: Better name?
struct AssembleRhsMap{A,B} <: Gridap.Fields.Map 
    agg_cells_local_dof_ids::A 
    agg_cells_rhs_contribs::B
end 

function Gridap.Fields.return_cache(m::AssembleRhsMap,aggregate_local_cells)
  cache_agg_cells_local_dof_ids=array_cache(m.agg_cells_local_dof_ids)
  cache_unassembled_rhs=array_cache(m.agg_cells_rhs_contribs)
  evaluate_result=Gridap.Arrays.CachedArray(eltype(eltype(m.agg_cells_rhs_contribs)),2)
  cache_agg_cells_local_dof_ids,cache_unassembled_rhs,evaluate_result
end

function Gridap.Fields.evaluate!(cache,m::AssembleRhsMap,aggregate_local_cells)
  cache_agg_cells_local_dof_ids,cache_unassembled_rhs,result=cache
  contrib = getindex!(cache_unassembled_rhs,m.agg_cells_rhs_contribs,1)

  max_local_dof_id=-1
  for (i,cell) in enumerate(aggregate_local_cells)
    current_cell_local_dof_ids = getindex!(cache_agg_cells_local_dof_ids,m.agg_cells_local_dof_ids,cell)
    for local_dof in current_cell_local_dof_ids
        max_local_dof_id=max(max_local_dof_id,local_dof)
    end
  end

  Gridap.Arrays.setsize!(result,(size(contrib,1),max_local_dof_id))
  result.array .= 0.0
  for (i,cell) in enumerate(aggregate_local_cells)
        current_cell_local_dof_ids = getindex!(cache_agg_cells_local_dof_ids,m.agg_cells_local_dof_ids,cell)
        contrib = getindex!(cache_unassembled_rhs,m.agg_cells_rhs_contribs,cell)
        for (j,local_dof) in enumerate(current_cell_local_dof_ids)
           result.array[:,local_dof] += contrib[:,j]
        end
  end
  result.array
end

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

agg_cells_local_dof_ids=compute_agg_cells_local_dof_ids(Ωagg_cell_dof_ids, aggregate_to_local_cells)
ass_rhs_map=AssembleRhsMap(agg_cells_local_dof_ids,agg_cells_rhs_contribs)
rhs=lazy_map(ass_rhs_map,aggregate_to_local_cells)

# TO-DO: optimize using our own optimized version Gridap.Fields.Map 
# of backslash that re-uses storage for lu factors among cells, etc.
dv_l2_proj_bb_dofs=lazy_map(\,lhs,rhs)

# Generate bb-wise array of fields. For each aggregate's bounding box,
# it provides the l2 projection of all basis functions in Ustd 
# restricted to the cells included in the bounding box of the aggregate   
dv_l2_proj_bb_array=lazy_map(Gridap.Fields.linear_combination,
                             dv_l2_proj_bb_dofs,
                             Gridap.CellData.get_data(vbb))

# Change domain of dv_l2_proj_bb_array from bb to agg_cells
dv_l2_proj_bb_array_agg_cells=lazy_map(Broadcasting(∘),
                                       lazy_map(Reindex(dv_l2_proj_bb_array),agg_cells_to_aggregate),
                                       ref_agg_cell_to_ref_bb_map)
dv_l2_proj_agg_cells=Gridap.CellData.GenericCellField(dv_l2_proj_bb_array_agg_cells,
                                                      Ωagg_cells,
                                                      ReferenceDomain()) 
du_l2_proj_agg_cells=Gridap.CellData.GenericCellField(lazy_map(transpose,dv_l2_proj_bb_array_agg_cells),
                                                      Ωagg_cells,
                                                      ReferenceDomain()) 

# Compute and assemble the bulk penalty stabilization term
# ∫( (dv-dv_l2_proj_agg_cells)*(du-du_l2_proj_agg_cells))*dΩ_agg_cells

γ = 10.0 # Interior bulk-penalty stabilization parameter
         # (@amartinhuertas no idea what a reasonable value is)

Ωbb  = Triangulation(aggregates_bounding_box_model)
dΩbb = Measure(Ωbb, degree)
h_U_array = get_array(∫(1.0)dΩbb)
h_U_array = lazy_map(Reindex(h_U_array), agg_cells_to_aggregate)
h_U = CellField(h_U_array, Ωagg_cells)

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

aggregate_dof_ids=compute_aggregate_dof_ids(Ωagg_cell_dof_ids,aggregate_to_cells)
agg_cells_to_aggregate_dof_ids=lazy_map(Reindex(aggregate_dof_ids),agg_cells_to_aggregate)

# Manually set up the arrays that collect_cell_matrix would return automatically
w = []
r = []
c = []

dv=get_fe_basis(Vstd)
du=get_trial_fe_basis(Ustd)

dv_du_mat_contribs=get_array(∫(γ*(1.0/h_U*h_U)*dv*du)*dΩagg_cells)
push!(w, dv_du_mat_contribs)
push!(r, Ωagg_cell_dof_ids)
push!(c, Ωagg_cell_dof_ids)

dv_proj_du_mat_contribs=get_array(∫(γ*(1.0/h_U*h_U)*(-1.0)*dv_l2_proj_agg_cells*du)*dΩagg_cells)
push!(w, dv_proj_du_mat_contribs)
push!(r, agg_cells_to_aggregate_dof_ids)
push!(c, Ωagg_cell_dof_ids)

proj_dv_du_mat_contribs=get_array(∫(γ*(1.0/h_U*h_U)*(-1.0)*dv*(du_l2_proj_agg_cells))*dΩagg_cells)
push!(w, proj_dv_du_mat_contribs)
push!(r, Ωagg_cell_dof_ids)
push!(c, agg_cells_to_aggregate_dof_ids)

proj_dv_proj_du_mat_contribs=get_array(∫(γ*(1.0/h_U*h_U)*dv_l2_proj_agg_cells*du_l2_proj_agg_cells)dΩagg_cells)
push!(w, proj_dv_proj_du_mat_contribs)
push!(r, agg_cells_to_aggregate_dof_ids)
push!(c, agg_cells_to_aggregate_dof_ids)

# Set up global projection matrix 
Ωcut = Triangulation(cutdisk,PHYSICAL)
dΩcut = Measure(Ωcut,degree)
a(u,v)=∫(v*u)dΩcut 
wrc=Gridap.FESpaces.collect_cell_matrix(Ustd,Vstd,a(du,dv))

assem=SparseMatrixAssembler(Ustd,Vstd)
Anostab=assemble_matrix(assem, wrc)
cond(Array(Anostab))

# Add the bulk penalty stabilization term to wrc
push!(wrc[1], w...)
push!(wrc[2], r...)
push!(wrc[3], c...)

Awithstab=assemble_matrix(assem, wrc)
cond(Array(Awithstab))

# Set up rhs global projection 
l(v)=∫(v*uex)dΩcut
b = assemble_vector(l, Vstd)

global_l2_proj_dofs = Awithstab\b
uh = FEFunction(Ustd, global_l2_proj_dofs)
eh = uex-uh
@assert sum(∫(eh*eh)*dΩcut) < 1.0e-12
