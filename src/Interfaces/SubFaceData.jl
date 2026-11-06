"""
    struct SubFaceData{Df,Dr,Dp,Tp,Tr,N} <: GridapType

Collection of `Df`-dimensional simplices (sub-faces) attached to the `Dr`-dimensional faces
of a background mesh, embedded in a `Dp`-dimensional space (`Df ≤ Dr ≤ Dp`):

- `cell_to_points`: vertices of each simplex (`Df+1` points per simplex);
- `cell_to_bgcell`: the background `Dr`-face containing each simplex;
- `point_to_coords`: physical coordinates of the points;
- `point_to_rcoords`: reference coordinates of the points in their background `Dr`-face;
- `cell_to_normal`: `nothing`, or the unit normal of each simplex (only if `Df == Dp-1`).

## Aliases: 

`SubCellData{Dr,Dp,Tp,Tr}` and `SubFacetData{Dp,T,Tn}` are now aliases.

    const SubCellData{Dr,Dp,Tp,Tr} = SubFaceData{Dr,Dr,Dp,Tp,Tr,Nothing}
    const SubFacetData{Dp,T,Tn} = SubFaceData{Df,Dp,Dp,T,T,Vector{VectorValue{Dp,Tn}}} where Df

## Constructors: 

    SubFaceData{Df}(cell_to_points, cell_to_bgcell, point_to_coords, point_to_rcoords[, cell_to_normal])
    SubCellData(cell_to_points, cell_to_bgcell, point_to_coords, point_to_rcoords)
    SubFacetData(facet_to_points, facet_to_normal, facet_to_bgcell, point_to_coords, point_to_rcoords)

"""
struct SubFaceData{Df,Dr,Dp,Tp,Tr,N} <: GridapType
  cell_to_points::Table{Int32,Vector{Int32},Vector{Int32}}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,Tp}}
  point_to_rcoords::Vector{Point{Dr,Tr}}
  cell_to_normal::N
  function SubFaceData{Df,Dr,Dp,Tp,Tr,N}(
    cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords,cell_to_normal
  ) where {Df,Dr,Dp,Tp,Tr,N}
    @assert 0 <= Df <= Dr <= Dp "SubFaceData{Df,Dr,Dp} requires 0 ≤ Df ≤ Dr ≤ Dp"
    @assert (N === Nothing) || (N <: Vector{<:VectorValue{Dp}} && Df == Dp-1) """
    The normals of a SubFaceData{Df,Dr,Dp} are either `nothing` or a Vector{VectorValue{Dp}},
    and can only be stored when Df == Dp-1.
    """
    new{Df,Dr,Dp,Tp,Tr,N}(
      cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords,cell_to_normal)
  end
end

const SubCellData{Dr,Dp,Tp,Tr} = SubFaceData{Dr,Dr,Dp,Tp,Tr,Nothing}
const SubFacetData{Dp,T,Tn} = SubFaceData{Df,Dp,Dp,T,T,Vector{VectorValue{Dp,Tn}}} where Df

function SubFaceData{Df,Dr,Dp,Tp,Tr,Nothing}(
  cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords) where {Df,Dr,Dp,Tp,Tr}
  SubFaceData{Df,Dr,Dp,Tp,Tr,Nothing}(
    cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords,nothing)
end

function SubFaceData{Df}(
  cell_to_points,
  cell_to_bgcell,
  point_to_coords::AbstractVector{<:Point{Dp,Tp}},
  point_to_rcoords::AbstractVector{<:Point{Dr,Tr}},
  cell_to_normal=nothing) where {Df,Dp,Tp,Dr,Tr}
  normals = _convert_normals(cell_to_normal)
  N = typeof(normals)
  SubFaceData{Df,Dr,Dp,Tp,Tr,N}(
    cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords,normals)
end

_convert_normals(::Nothing) = nothing
_convert_normals(n::AbstractVector{<:VectorValue}) = convert(Vector{eltype(n)},n)

function SubCellData(
  cell_to_points,
  cell_to_bgcell,
  point_to_coords::AbstractVector{<:Point},
  point_to_rcoords::AbstractVector{<:Point{Dr}}) where Dr
  SubFaceData{Dr}(cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords)
end

function SubFacetData(
  facet_to_points,
  facet_to_normal::AbstractVector{<:VectorValue},
  facet_to_bgcell,
  point_to_coords::AbstractVector{<:Point{Dp}},
  point_to_rcoords) where Dp
  SubFaceData{Dp-1}(
    facet_to_points,facet_to_bgcell,point_to_coords,point_to_rcoords,facet_to_normal)
end

# Former field names of SubFacetData (deprecated, kept for downstream packages)

function Base.getproperty(st::SubFaceData,sym::Symbol)
  if sym === :facet_to_points
    getfield(st,:cell_to_points)
  elseif sym === :facet_to_bgcell
    getfield(st,:cell_to_bgcell)
  elseif sym === :facet_to_normal
    getfield(st,:cell_to_normal)
  else
    getfield(st,sym)
  end
end

Base.propertynames(st::SubFaceData,private::Bool=false) = fieldnames(typeof(st))

Gridap.Geometry.num_cell_dims(::SubFaceData{Df}) where Df = Df
Gridap.Geometry.num_cell_dims(::Type{<:SubFaceData{Df}}) where Df = Df
Gridap.Geometry.num_point_dims(::SubFaceData{Df,Dr,Dp}) where {Df,Dr,Dp} = Dp
Gridap.Geometry.num_point_dims(::Type{<:SubFaceData{Df,Dr,Dp}}) where {Df,Dr,Dp} = Dp

# Restriction

"""
    SubFaceData(st::SubFaceData, newcells)
    SubFaceData(st::SubFaceData, newcells, orientation)

Restriction of `st` to the simplices `newcells` (all the points are kept). 
With `orientation` (`±1` per new simplex), the stored normals are 
multiplied by it (requires stored normals).
"""
function SubFaceData(st::SubFaceData{Df},newcells::AbstractVector{<:Integer}) where Df
  cell_to_points = Table(lazy_map(Reindex(st.cell_to_points),newcells))
  cell_to_bgcell = st.cell_to_bgcell[newcells]
  if isnothing(st.cell_to_normal)
    cell_to_normal = nothing
  else
    cell_to_normal = st.cell_to_normal[newcells]
  end
  return SubFaceData{Df}(
    cell_to_points, cell_to_bgcell,
    st.point_to_coords, st.point_to_rcoords, cell_to_normal
  )
end

function SubFaceData(
  st::SubFaceData{Df,Dr,Dp},
  newcells::AbstractVector{<:Integer},
  orientation::AbstractVector) where {Df,Dr,Dp}
  @notimplementedif isnothing(st.cell_to_normal) "Orientation requires stored normals"
  cell_to_points = Table(lazy_map(Reindex(st.cell_to_points),newcells))
  cell_to_bgcell = st.cell_to_bgcell[newcells]
  cell_to_normal = st.cell_to_normal[newcells] .* orientation
  SubFaceData{Df}(
    cell_to_points, cell_to_bgcell,
    st.point_to_coords, st.point_to_rcoords, cell_to_normal
  )
end

SubCellData(st::SubCellData,newcells::AbstractVector{<:Integer}) = SubFaceData(st,newcells)

function SubFacetData(
  st::SubFacetData,newfacets::AbstractVector{<:Integer},orientation::AbstractVector)
  return SubFaceData(st,newfacets,orientation)
end

# In-place concatenation

function Base.empty(st::SubFaceData{Df,Dr,Dp,Tp,Tr}) where {Df,Dr,Dp,Tp,Tr}
  normals = isnothing(st.cell_to_normal) ? nothing : empty(st.cell_to_normal)
  return SubFaceData{Df}(
    Table(Int32[],Int32[1,]), Int32[], Point{Dp,Tp}[], Point{Dr,Tr}[], normals
  )
end

function Base.append!(a::SubFaceData,b::SubFaceData)
  @assert num_cell_dims(a) == num_cell_dims(b)
  @assert isnothing(a.cell_to_normal) == isnothing(b.cell_to_normal)

  o = length(a.point_to_coords)

  if !isnothing(a.cell_to_normal)
    append!(a.cell_to_normal, b.cell_to_normal)
  end
  append!(a.cell_to_bgcell, b.cell_to_bgcell)
  append!(a.point_to_coords, b.point_to_coords)
  append!(a.point_to_rcoords, b.point_to_rcoords)

  nini = length(a.cell_to_points.data)+1
  append!(a.cell_to_points.data, b.cell_to_points.data)
  nend = length(a.cell_to_points.data)
  append_ptrs!(a.cell_to_points.ptrs,b.cell_to_points.ptrs)
  for i in nini:nend
    a.cell_to_points.data[i] += o
  end

  return a
end

function merge_sub_face_data(ls_to_subfacets,ls_to_ls_to_facet_to_inout)
  fst = empty(first(ls_to_subfacets))
  ls_to_facet_to_inout = [ empty(first(ls_to_f_to_i))  for ls_to_f_to_i in ls_to_ls_to_facet_to_inout  ]
  for (i,fst_i) in enumerate(ls_to_subfacets)
    append!(fst,fst_i)
    for (j,facet_to_inout) in enumerate(ls_to_ls_to_facet_to_inout[i])
      append!(ls_to_facet_to_inout[j],facet_to_inout)
    end
  end
  return fst, ls_to_facet_to_inout
end

# Grid, reference map of the glue, visualization

# Stored normals become the grid's facet normals (`get_facet_normal`)
function UnstructuredGrid(st::SubFaceData{Df}) where Df
  reffe = LagrangianRefFE(Float64,Simplex(Val{Df}()),1)
  cell_types = fill(Int8(1),length(st.cell_to_points))
  if isnothing(st.cell_to_normal)
    normals = nothing
  else
    normals = lazy_map(constant_field, st.cell_to_normal)
  end
  return UnstructuredGrid(
    st.point_to_coords, st.cell_to_points, [reffe,], 
    cell_types, NonOriented(), normals
  )
end

# Map from the reference simplex to the reference space of the background Dr-face
function _setup_cell_ref_map(st::SubFaceData,grid)
  cell_to_points = st.cell_to_points
  point_to_rcoords = st.point_to_rcoords
  cell_to_rcoords = lazy_map(Broadcasting(Reindex(point_to_rcoords)),cell_to_points)
  ctype_to_reffe = get_reffes(grid)
  cell_to_ctype = get_cell_type(grid)
  @notimplementedif length(ctype_to_reffe) != 1
  reffe = first(ctype_to_reffe)
  cell_to_shapefuns = expand_cell_data([get_shapefuns(reffe)],cell_to_ctype)
  cell_to_ref_map = lazy_map(linear_combination,cell_to_rcoords,cell_to_shapefuns)
  return cell_to_ref_map
end

function Visualization.visualization_data(
  st::SubFaceData{Df,Dr},filename::String;celldata=Dict()) where {Df,Dr}
  ug = UnstructuredGrid(st)
  degree = 0
  trian = GenericTriangulation(ug)
  quad = CellQuadrature(trian,degree)
  measure = integrate(1,quad)
  newcelldata = Pair{String,Any}["bgcell"=>st.cell_to_bgcell, (Df == Dr ? "dV" : "dS")=>measure]
  if !isnothing(st.cell_to_normal)
    push!(newcelldata,"normal"=>st.cell_to_normal)
  end
  _celldata = Dict()
  for (k,v) in celldata
    _celldata[k] = v
  end
  for (k,v) in newcelldata
    _celldata[k] = v
  end
  return (Visualization.VisualizationData(ug,filename,celldata=_celldata),)
end
