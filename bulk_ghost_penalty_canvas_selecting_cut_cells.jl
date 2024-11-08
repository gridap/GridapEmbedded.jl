using Gridap
using GridapEmbedded
using FillArrays
using LinearAlgebra

include("BulkGhostPenaltyAssembleMaps.jl")

# Manufactured solution
order = 0
uex(x) = -VectorValue(2*x[1],2*x[2])
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

#INVESTIGATING: Number of aggregates, cells and bg cells.
@assert length(aggregates) == n*n
num_cell = length(aggregates)
using StatsBase
counts = countmap(aggregates,alg=:dict)
agg_counts = filter(((k,v),) -> v > 1 && k != 0, counts)
num_aggregates = length(agg_counts)
num_cells_in_aggregates = 0
for (k,v) in agg_counts
    println(v)
    num_cells_in_aggregates +=v
end
num_cells_not_in_aggregates = counts[0]
non_agg_counts = filter(((k,v),) -> v <= 1 || k==0, counts)
num_cells_not_in_aggregates = 0
for (k,v) in non_agg_counts
    println(v)
    num_cells_not_in_aggregates +=v
end
num_ext_cells = counts[0]
@assert num_cells_in_aggregates + num_cells_not_in_aggregates == num_cell
println("Out of the $num_cell cells, $num_cells_in_aggregates belong to aggregates and $num_cells_not_in_aggregates are not part of aggregates. Note that we have $num_ext_cells exterior cells")
println("Assuming that there is one root per aggregate, we have $(num_cells_in_aggregates-num_aggregates) cut cells and $num_aggregates roots in the $num_aggregates aggregates cotntaining $num_cells_in_aggregates in total.")

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
  Creates an array of arrays with as many entries 
  as interior cells that are not part of any aggegrate. 
  For each interior cell, the array 
  contains the global cell IDs of that cells in the background 
  model that belong to the same interior cell

  TO-DO: with efficiency in mind we may want to store this 
         array of arrays as a Gridap.Arrays.Table.
"""
function setup_int_nonagg_cell_to_cells(aggregates)
  
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
  interior_cell_to_cells=Vector{Vector{Int}}()
  current_interior_cell=1
  for (i,agg) in enumerate(aggregates)
    if agg>0
        if (size_aggregates[agg]==1)
            if !haskey(touched,agg)
                push!(interior_cell_to_cells,[i])
                touched[agg]=current_interior_cell
                current_interior_cell+=1
            else 
                push!(interior_cell_to_cells[touched[agg]],i)
            end
        end 
    end
  end
  interior_cell_to_cells
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

int_nonagg_cell_to_cells = setup_int_nonagg_cell_to_cells(aggregates)
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
                                    Ω_agg_cells,
                                       agg_cells_to_aggregate)
    @assert num_cells(Ω_agg_cells)==length(ref_agg_cell_to_ref_bb_map)
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
    
    bb_basis_array_to_Ω_agg_cells_array = lazy_map(Reindex(bb_basis_array),agg_cells_to_aggregate)
    bb_basis_array_to_Ω_agg_cells_array = lazy_map(Broadcasting(∘),
                                                  bb_basis_array_to_Ω_agg_cells_array,
                                                  ref_agg_cell_to_ref_bb_map)
    if (bb_basis_style==Gridap.FESpaces.TrialBasis())
        # Add transpose 
        bb_basis_array_to_Ω_agg_cells_array=lazy_map(transpose, bb_basis_array_to_Ω_agg_cells_array)
    end    

    Gridap.CellData.GenericCellField(bb_basis_array_to_Ω_agg_cells_array,
                                    Ω_agg_cells,
                                     ReferenceDomain())
end

# Generate an array with the global IDs of the cells 
# that belong to an aggregrate. From now on, we will 
# use the terminology "agg_cells" to refer to those 
# cells of the background model that belong to an aggregate 
# (i.e., they can be either cut or interior cells)
agg_cells=Vector{Int}()
for cells in aggregate_to_cells
    append!(agg_cells,cells)
end

# Generate an array with the global IDs of the cells 
# that belong to the interior, yet do not belong to the 
# aggregate. We will use "int_nonagg_cells" to refer to those 
# cells of the background model that belong to the interior 
# but are not part of any of the aggregates.
int_nonagg_cells=Vector{Int}()
for cell in int_nonagg_cell_to_cells
    append!(int_nonagg_cells,cell)
end

#INSPECTION
# agg_cells
# int_nonagg_cells
# aggregates

## Triangulation
Ωbg  = Triangulation(bgmodel)
Ω    = Triangulation(cutgeo,PHYSICAL)
Ωact = Triangulation(cutgeo,ACTIVE)
Ωcut = Triangulation(cutgeo)

# Interior cells
Ωbg_agg_cells        = view(Ωbg,agg_cells)            # cells in aggregates
int_cells = Ω_agg_cells.parent.b.tface_to_mface 
# construct the interior aggregate cells, in other words the root cells
root_cells=Vector{Int}()
tester_int_nonagg_cells=Vector{Int}()
for cell in int_cells
    if cell ∈ int_nonagg_cells
        append!(tester_int_nonagg_cells,cell)
    else
        append!(root_cells,cell)
    end
end
@assert(tester_int_nonagg_cells==int_nonagg_cells)
# Cut cells
cut_cells=Vector{Int}()
tester_root_cells=Vector{Int}()
for cell in agg_cells
    if cell ∈ root_cells
        append!(tester_root_cells,cell)
    else
        append!(cut_cells,cell)
    end
end
@assert(tester_root_cells==root_cells)


# TODO Potentially useful?:
# Ω_agg_cells.parent.a.subgrid.cell_node_ids.data == Ω.a.subgrid.cell_node_ids.data # vector containing cell ids 

# Subdomains in background mesh (aggegrate cells are already defined)
Ωbg_int_cells        = view(Ωbg,int_cells)            # cells in interior
Ωbg_int_nonagg_cells = view(Ωbg,int_nonagg_cells)     # cells in interior but not in aggregate
Ωbg_root_cells       = view(Ωbg,root_cells)           # root cells: in interior and part of aggregate
Ωbg_cut_cells        = view(Ωbg,cut_cells)            # cut cells (part of aggregate)

# Subdomains in physical mesh
Ω_agg_cells        = view(Ω,agg_cells)            # cells in aggregates
Ω_int_cells        = view(Ω,int_cells)            # cells in interior
Ω_int_nonagg_cells = view(Ω,int_nonagg_cells)     # cells in interior but not in aggregate
Ω_root_cells       = view(Ω,root_cells)           # root cells: in interior and part of aggregate
Ω_cut_cells        = view(Ω,cut_cells)            # cut cells (part of aggregate)

## VISUALISATION:
writevtk(Ωbg,"trian_bg")
writevtk(Ωbg_agg_cells,"trian_bg_agg_cells")
writevtk(Ωbg_int_cells,"trian_bg_int_cells")
writevtk(Ωbg_int_nonagg_cells,"trian_bg_int_nonagg_cells")
writevtk(Ωbg_root_cells,"trian_bg_root_cells")
writevtk(Ωbg_cut_cells,"trian_bg_cut_cells")

writevtk(Ωact,"trian_act")
writevtk(Ω,"trian_phys")
#TODO: WHY DOESN'T THIS WORK?
# writevtk(Ω_agg_cells,"trian_phys_agg_cells")
# writevtk(Ω_int_cells,"trian_phys_int_cells")
# writevtk(Ω_int_nonagg_cells,"trian_phys_int_nonagg_cells")
# writevtk(Ω_root_cells,"trian_phys_root_cells")
# writevtk(Ω_cut_cells,"trian_phys_cut_coot_cells")
degree = 2*2*order
dΩ    = Measure(Ω,degree)
dΩact = Measure(Ωact,degree)
dΩbg_agg_cells  = Measure(Ωbg_agg_cells,degree)
dΩbg_cut_cells  = Measure(Ωbg_cut_cells,degree)
dΩbg_root_cells = Measure(Ωbg_root_cells,degree)


dΩ_agg_cells  = Measure(Ω_agg_cells,degree)
dΩcut_cells  = Measure(Ω_cut_cells,degree)
dΩroot_cells = Measure(Ω_root_cells,degree)

# TESTING
area_phys_domain = sum(∫(1.0)dΩ)
area_act_domain  = sum(∫(1.0)dΩact)
area_agg_cells   = sum(∫(1.0)dΩbg_agg_cells)
area_cut_cells   = sum(∫(1.0)dΩbg_cut_cells)
area_root_cells  = sum(∫(1.0)dΩbg_root_cells)
@assert(area_agg_cells ≈ area_cut_cells + area_root_cells)

# TODO: not possible to integrate over the physical part of the domain, e.g. 
area_phys_agg_cells  = sum(∫(1.0)dΩ_agg_cells)
area_phys_cut_cells  = sum(∫(1.0)dΩ_cut_cells)
area_phys_root_cells = sum(∫(1.0)dΩ_root_cells)