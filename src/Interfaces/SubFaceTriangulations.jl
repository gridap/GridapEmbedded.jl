# Implementation of the Gridap.Triangulation interface

"""
    struct SubFaceTriangulation{Df,Dp,T,A,S<:SubFaceData,G} <: Triangulation{Df,Dp}

Triangulation of a [`SubFaceData`](@ref), glued to the background `Dr`-faces
(`get_glue(trian,Val(Dr))`). The stored normals of the data (if any) are the facet normals
of its grid (`get_facet_normal`); without them `get_facet_normal` errors. Aliases (same
parameters as the former structs):

- `SubCellTriangulation{Dc,Dp,T,A}`: on a `SubCellData`;
- `SubFacetTriangulation{Dc,Dp,T,A}`: on a `SubFacetData`.
"""
struct SubFaceTriangulation{Df,Dp,T,A,S<:SubFaceData,G<:UnstructuredGrid{Df,Dp,T}} <: Triangulation{Df,Dp}
  subfaces::S
  bgmodel::A
  subgrid::G
  function SubFaceTriangulation(
    subfaces::SubFaceData{Df,Dr,Dp,T},bgmodel::DiscreteModel) where {Df,Dr,Dp,T}
    subgrid = UnstructuredGrid(subfaces)
    A = typeof(bgmodel)
    S = typeof(subfaces)
    G = typeof(subgrid)
    new{Df,Dp,T,A,S,G}(subfaces,bgmodel,subgrid)
  end
end

"""
    const SubCellTriangulation{Dc,Dp,T,A} = SubFaceTriangulation{Dc,Dp,T,A,<:SubCellData{Dc,Dp,T}}

A triangulation for subcells.
"""
const SubCellTriangulation{Dc,Dp,T,A} = SubFaceTriangulation{Dc,Dp,T,A,<:SubCellData{Dc,Dp,T}}

"""
    const SubFacetTriangulation{Dc,Dp,T,A} = SubFaceTriangulation{Dc,Dp,T,A,<:SubFacetData{Dp,T}}

A triangulation for subfacets (pieces of the interface inside background cells).
"""
const SubFacetTriangulation{Dc,Dp,T,A} = SubFaceTriangulation{Dc,Dp,T,A,<:SubFacetData{Dp,T}}

SubCellTriangulation(st::SubCellData,bgmodel::DiscreteModel) = SubFaceTriangulation(st,bgmodel)
SubFacetTriangulation(st::SubFacetData,bgmodel::DiscreteModel) = SubFaceTriangulation(st,bgmodel)

# Former field names (deprecated, kept for downstream packages)
function Base.getproperty(a::SubFaceTriangulation,sym::Symbol)
  if sym === :subcells || sym === :subfacets
    getfield(a,:subfaces)
  else
    getfield(a,sym)
  end
end

Base.propertynames(a::SubFaceTriangulation,private::Bool=false) = fieldnames(typeof(a))

function get_background_model(a::SubFaceTriangulation)
  a.bgmodel
end

function get_grid(a::SubFaceTriangulation)
  a.subgrid
end

function get_glue(a::SubFaceTriangulation,::Val{D}) where D
  Dr = _ref_dim(a.subfaces)
  if D != Dr
    return nothing
  end
  tface_to_mface = a.subfaces.cell_to_bgcell
  tface_to_mface_map = _setup_cell_ref_map(a.subfaces,a.subgrid)
  FaceToFaceGlue(tface_to_mface,tface_to_mface_map,nothing)
end

_ref_dim(::SubFaceData{Df,Dr}) where {Df,Dr} = Dr

function move_contributions(scell_to_val::AbstractArray,strian::SubFaceTriangulation)
  model = get_background_model(strian)
  ncells = num_cells(model)
  cell_to_touched = fill(false,ncells)
  scell_to_cell = strian.subfaces.cell_to_bgcell
  cell_to_touched[scell_to_cell] .= true
  Ωa = Triangulation(model,cell_to_touched)
  acell_to_val = move_contributions(scell_to_val,strian,Ωa)
  acell_to_val, Ωa
end

# Compute the active model
# To do this, we need to glue together the sub-faces, which unfortunately requires comparing
# the point coordinates...
function compute_active_model(trian::SubFaceTriangulation)
  subgrid = trian.subgrid
  subfaces = trian.subfaces
  cell_to_uids, uid_to_point = consistent_facet_to_points(
    subfaces.cell_to_points, subfaces.point_to_coords
  )
  topo = UnstructuredGridTopology(
    subgrid, cell_to_uids, uid_to_point
  )
  return UnstructuredDiscreteModel(subgrid,topo,FaceLabeling(topo))
end

function consistent_facet_to_points(
  facet_to_points::Table, point_to_coords::Vector
)
  f(pt::VectorValue) = VectorValue(round.(pt.data;sigdigits=12))
  f(id::Integer) = f(point_to_coords[id])

  # Create a list of the unique points composing the facets
  npts = length(point_to_coords)
  nfaces = length(facet_to_points)
  touched = zeros(Bool,npts)
  for face in 1:nfaces
    pts = view(facet_to_points,face)
    touched[pts] .= true
  end
  touched_ids = findall(touched)
  unique_ids = unique(f,touched_ids)

  # Create a mapping from the old point ids to the new ones
  touched_to_uid = collect(Int32,indexin(f.(touched_ids),f.(unique_ids)))
  point_to_uid = extend(touched_to_uid,PosNegPartition(touched_ids,npts))

  facet_to_uids = Table(
    collect(Int32,lazy_map(Reindex(point_to_uid),facet_to_points.data)),
    facet_to_points.ptrs
  )
  return facet_to_uids, unique_ids
end
