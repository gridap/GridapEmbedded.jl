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

"""
  Creates an array containing the cell IDs of a selection of the background cells. The cell type filter IN (-1) selects the interior cells, OUT (1) the exterior cells, and GridapEmbedded.Interfaces.CUT (0) the cut cells.
  TO-DO: publish CUT, so that GridapEmbedded.Interfaces.CUT can be shortened to CUT?
  TO-DO: be careful with using restrict_cells(cutgeo,GridapEmbedded.Interfaces.CUT) to replace flatten(restrict_aggregate_to_cells(cutgeo,aggregate_to_cells,GridapEmbedded.Interfaces.CUT))

"""
function restrict_cells(cutgeo,cell_type_filter)
  restricted_cells=Int64[]
  cell_to_inoutcut=compute_bgcell_to_inoutcut(cutgeo,cutgeo.geo)
  for cell in 1:length(cell_to_inoutcut)
    if cell_to_inoutcut[cell] == cell_type_filter
        push!(restricted_cells, cell)
    end
  end
  restricted_cells
end

"""
  Creates an array of arrays with as many entries 
  as aggregates. For each aggregate, the array 
  contains a selection of the global cell IDs of 
  the cells in the background model that belong to 
  the same aggregate. The cell type filter IN (-1)
  selects the interior cells, OUT (1) the exterior
  cells, and GridapEmbedded.Interfaces.CUT (0) the
  cut cells.
  TO-DO: publish CUT, so that GridapEmbedded.Interfaces.CUT can be shortened to CUT?
"""
function restrict_aggregate_to_cells(cutgeo,aggregate_to_cells,cell_type_filter)
  aggregate_to_restricted_cells=Vector{Int}[]
  cell_to_inoutcut=compute_bgcell_to_inoutcut(cutgeo,cutgeo.geo)
  for agg in aggregate_to_cells
    restricted_cells = Int64[]
    for cell in agg
        if cell_to_inoutcut[cell] == cell_type_filter
            push!(restricted_cells, cell)
        end
    end
    push!(aggregate_to_restricted_cells,restricted_cells)
  end
  aggregate_to_restricted_cells
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

"""
  Flattens an array of arrays
"""
function flatten(array_of_arrays)
    vcat(array_of_arrays...)
 end

"""
  Generate an array that given the local ID of an "agg_cell"
  returns the ID of the aggregate to which it belongs
  (i.e., flattened version of aggregate_to_cells)
    
  
  TO-DO: perhaps merge this function with , 
  setup_cut_cells_in_agg_cells_to_aggregate

"""
function setup_cells_to_aggregate(aggregate_to_cells)
    cells_to_aggregate=Vector{Int}()
    for (i,cells) in enumerate(aggregate_to_cells)
        for _ in cells 
            push!(cells_to_aggregate,i)
        end
    end
    cells_to_aggregate 
end

"""
  Renumbers the global cell IDs to local cell IDs in the array of arrays
  with as many entries as aggregates. For each aggregate, the array 
  now contains the local cell IDs of that cells in the background 
  model that belong to the same aggregate

  TO-DO: with efficiency in mind we may want to store this 
         array of arrays as a Gridap.Arrays.Table.
"""
function setup_aggregate_to_local_cells(aggregate_to_cells)
    aggregate_to_local_cells=deepcopy(aggregate_to_cells)
    current_local_cell=1
    for (i,cells) in enumerate(aggregate_to_local_cells)
        for j in 1:length(cells)
            cells[j]=current_local_cell
            current_local_cell+=1
        end
    end
    aggregate_to_local_cells
end

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

"""
  Define mapping ref_agg_cell_to_ref_bb_map: K_ref -> K -> bb -> bb_ref
"""
function setup_ref_agg_cell_to_ref_bb_map(aggregates_bounding_box_model,agg_cells_to_aggregate,ref_agg_cell_to_agg_cell_map)
    bb_to_ref_bb=lazy_map(Gridap.Fields.inverse_map,get_cell_map(aggregates_bounding_box_model))
    bb_to_ref_bb_agg_cells=lazy_map(Reindex(bb_to_ref_bb),agg_cells_to_aggregate)  
    ref_agg_cell_to_ref_bb_map=
      lazy_map(Broadcasting(∘),bb_to_ref_bb_agg_cells,ref_agg_cell_to_agg_cell_map)
end

"""
  Compute local dof ids of the aggregated cells
"""
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

"""
  Returns the dof ids for the aggregates
"""
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