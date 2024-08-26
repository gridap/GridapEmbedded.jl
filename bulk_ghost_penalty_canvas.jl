using Gridap
using GridapEmbedded
using FillArrays

# Manufactured solution
order = 1
u(x) = x[1] + x[2]^order
∇u(x) = ∇(u)(x)
f(x) = -Δ(u)(x)
ud(x) = u(x)

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
  touched=Dict{Int,Int}()
  aggregate_to_cells=Vector{Vector{Int}}()
  current_aggregate=1
  for (i,agg) in enumerate(aggregates)
    if agg>0
        if !haskey(touched,agg)
            push!(aggregate_to_cells,[i])
            touched[agg]=current_aggregate
            current_aggregate+=1
        else 
            push!(aggregate_to_cells[touched[agg]],i)
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
order = 1
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
# (i.e., flattened version of aggregrate_to_cells)
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
degree=2*order
dΩagg_cells = Measure(Ωagg_cells,degree)
order=1
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

function Gridap.Fields.return_cache(m::AssembleLhsMap,cells)
    cache_unassembled_lhs=array_cache(m.agg_cells_lhs_contribs)
    evaluate_result=Gridap.Arrays.CachedArray(eltype(eltype(m.agg_cells_lhs_contribs)),2)
    cache_unassembled_lhs,evaluate_result
end
  
function Gridap.Fields.evaluate!(cache,m::AssembleLhsMap,cells)
    cache_unassembled_lhs,result=cache
    contrib = getindex!(cache_unassembled_lhs,m.agg_cells_lhs_contribs,1)
    Gridap.Arrays.setsize!(result,(size(contrib,1),size(contrib,2)))
    result.array .= 0.0
    for (i,cell) in enumerate(cells)
        contrib = getindex!(cache_unassembled_lhs,m.agg_cells_lhs_contribs,cell)
        result.array .= contrib
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
    agg_cells_dof_ids::A 
    agg_cells_rhs_contribs::B
end 

function Gridap.Fields.return_cache(m::AssembleRhsMap,cells)
  cache_agg_cells_dof_ids=array_cache(m.agg_cells_dof_ids)
  cache_unassembled_rhs=array_cache(m.agg_cells_rhs_contribs)
  evaluate_result=Gridap.Arrays.CachedArray(eltype(eltype(m.agg_cells_rhs_contribs)),2)
  cache_agg_cells_dof_ids,cache_unassembled_rhs,evaluate_result
end

function Gridap.Fields.evaluate!(cache,m::AssembleRhsMap,cells)
  cache_cell_dof_ids,cache_unassembled_rhs,result=cache

  # TO-DO: this is a naive/not efficient implementation to 
  #        compute a local numbering of dofs of those cells 
  #        which belong to current aggregate. We should 
  #        pre-compute this numbering outside and store it in an 
  #        array of arrays (i.e., Gridap.Arrays.Table) 
  global_to_local_dofs = Dict{Int32,Int32}()
  current=1
  for (i,cell) in enumerate(cells)
        current_cell_dof_ids = getindex!(cache_cell_dof_ids,
                                         m.agg_cells_dof_ids,
                                         cell) 
        for dof in current_cell_dof_ids
          if !haskey(global_to_local_dofs,dof)
            global_to_local_dofs[dof]=current
            current+=1
          end 
        end 
  end

  contrib = getindex!(cache_unassembled_rhs,m.agg_cells_rhs_contribs,1)
  Gridap.Arrays.setsize!(result,(size(contrib,1),current-1))
  result.array .= 0.0
  for (i,cell) in enumerate(cells)
        current_cell_dof_ids = getindex!(cache_cell_dof_ids,m.agg_cells_dof_ids,cell)
        contrib = getindex!(cache_unassembled_rhs,m.agg_cells_rhs_contribs,cell)
        for (j,dof) in enumerate(current_cell_dof_ids)
           result.array[:,global_to_local_dofs[dof]] += contrib[:,j]
        end
  end
  result.array
end

ass_rhs_map=AssembleRhsMap(Ωagg_cell_dof_ids,agg_cells_rhs_contribs)
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

# TO-DO: how to implement? ...
# ∫( (dv-dv_l2_proj_agg_cells)*(du-du_l2_proj_agg_cells))*dΩ_agg_cells

dv=get_fe_basis(Vstd)
du=get_trial_fe_basis(Ustd)

get_array(∫(dv*du)*dΩagg_cells)[1] #- 
  get_array(∫(dv*du_l2_proj_agg_cells)*dΩagg_cells)[1] #- 
    get_array(∫(dv_l2_proj_agg_cells*du)*dΩagg_cells)[1] #+ 
       get_array(∫(dv_l2_proj_agg_cells*du_l2_proj_agg_cells)*dΩagg_cells)[1]
