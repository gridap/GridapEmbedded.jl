"""
    struct ImmersedFaceData{Df,Dp,T}

The cutfaces of the `Df`-ifaces: contains a `SubFaceData{Df,Dp,Dp}`, plus 
`point_to_ircoords` (reference coordinates of the iface the cutface lies on).
"""
struct ImmersedFaceData{Df,Dp,T} <: GridapType
  subfaces::SubFaceData{Df,Dp,Dp,T,T,Nothing}
  point_to_ircoords::Vector{Point{Df,T}}
end

function Base.getproperty(st::ImmersedFaceData,sym::Symbol)
  if sym in (:subfaces,:point_to_ircoords)
    getfield(st,sym)
  else
    getproperty(getfield(st,:subfaces),sym)
  end
end

Base.propertynames(st::ImmersedFaceData,private::Bool=false) =
  (:subfaces,:point_to_ircoords,fieldnames(typeof(st.subfaces))...)

Gridap.Geometry.num_cell_dims(::ImmersedFaceData{Df}) where Df = Df
Gridap.Geometry.num_point_dims(::ImmersedFaceData{Df,Dp}) where {Df,Dp} = Dp

"""
    ImmersedFaceData(st::ImmersedFaceData, newcells)

The cutfaces `newcells` of `st` and only their points (renumbered): every point belongs to a
cutface.
"""
function ImmersedFaceData(st::ImmersedFaceData{Df,Dp,T},newcells::AbstractVector{<:Integer}) where {Df,Dp,T}
  cell_to_oldpoints = Table(lazy_map(Reindex(st.cell_to_points),newcells))
  newpoint_to_oldpoint = unique(cell_to_oldpoints.data)
  oldpoint_to_newpoint = zeros(Int32,length(st.point_to_coords))
  oldpoint_to_newpoint[newpoint_to_oldpoint] .= 1:length(newpoint_to_oldpoint)
  cell_to_points = Table(oldpoint_to_newpoint[cell_to_oldpoints.data],cell_to_oldpoints.ptrs)
  subfaces = SubFaceData{Df}(
    cell_to_points,st.cell_to_bgcell[newcells],
    st.point_to_coords[newpoint_to_oldpoint],st.point_to_rcoords[newpoint_to_oldpoint])
  return ImmersedFaceData(subfaces,st.point_to_ircoords[newpoint_to_oldpoint])
end

# The D-cutfaces as subcells, the (D-1)-cutfaces as subfacets (with a normal per cutface):
# the iface reference coordinates are dropped.
SubCellData(st::ImmersedFaceData{D,D}) where D = st.subfaces
function SubFacetData(st::ImmersedFaceData{Df,D},cell_to_normal::AbstractVector) where {Df,D}
  s = st.subfaces
  SubFacetData(s.cell_to_points,collect(cell_to_normal),s.cell_to_bgcell,s.point_to_coords,s.point_to_rcoords)
end

UnstructuredGrid(st::ImmersedFaceData) = UnstructuredGrid(st.subfaces)

SubFaceTriangulation(st::ImmersedFaceData,bgmodel::DiscreteModel) = SubFaceTriangulation(st.subfaces,bgmodel)

"""
    struct ImmersedDiscreteModel{Dc,Dp} <: DiscreteModel{Dc,Dp}

Represents an immersed model `imodel` cut by the cells of `bgmodel`.
Created by `cut(cutter, bgmodel, imodel)`. 

Contains:

- `cutfaces`: tuple, entry d+1 the `ImmersedFaceData{d,Dp}` of the cutfaces of the d-ifaces,
  or `nothing` if dimension d was not cut;
- `iface_to_cutfaces`: entry d+1 the `Table` d-iface → cutfaces (undefined if not cut);
- `metadata`: how the cut was made.
"""
struct ImmersedDiscreteModel{Dc,Dp,A,B,C,M} <: DiscreteModel{Dc,Dp}
  imodel::A
  bgmodel::B
  cutfaces::C
  iface_to_cutfaces::Vector{Table{Int32,Vector{Int32},Vector{Int32}}}
  metadata::M
  function ImmersedDiscreteModel(
    imodel::DiscreteModel{Dc,Dp},bgmodel::DiscreteModel{Dc,Dp},cutfaces::Tuple,
    iface_to_cutfaces::Vector{Table{Int32,Vector{Int32},Vector{Int32}}},metadata) where {Dc,Dp}
    @check length(cutfaces) == Dc+1 && length(iface_to_cutfaces) == Dc+1
    A, B, C, M = typeof(imodel), typeof(bgmodel), typeof(cutfaces), typeof(metadata)
    new{Dc,Dp,A,B,C,M}(imodel,bgmodel,cutfaces,iface_to_cutfaces,metadata)
  end
end

get_grid(m::ImmersedDiscreteModel) = get_grid(m.imodel)
Gridap.Geometry.get_grid_topology(m::ImmersedDiscreteModel) = get_grid_topology(m.imodel)
Gridap.Geometry.get_face_labeling(m::ImmersedDiscreteModel) = get_face_labeling(m.imodel)

"""
    get_imodel(m::ImmersedDiscreteModel), get_bgmodel(m::ImmersedDiscreteModel)
"""
get_imodel(m::ImmersedDiscreteModel) = m.imodel
get_bgmodel(m::ImmersedDiscreteModel) = m.bgmodel

"""
    get_cutfaces(m::ImmersedDiscreteModel, d) -> ImmersedFaceData{d}

The cutfaces of the d-ifaces (error if dimension d was not cut). Type-stable with `Val(d)`.
"""
get_cutfaces(m::ImmersedDiscreteModel,d::Integer) = get_cutfaces(m,Val(d))
function get_cutfaces(m::ImmersedDiscreteModel,::Val{d}) where d
  data = m.cutfaces[d+1]
  isnothing(data) && error("Dimension $d of the ImmersedDiscreteModel was not cut")
  return data
end

"""
    get_iface_to_cutfaces(m::ImmersedDiscreteModel, d) -> Table

d-iface → its cutfaces (ids in `get_cutfaces(m,d)`).
"""
function get_iface_to_cutfaces(m::ImmersedDiscreteModel,d::Integer)
  isassigned(m.iface_to_cutfaces,d+1) || error("Dimension $d of the ImmersedDiscreteModel was not cut")
  return m.iface_to_cutfaces[d+1]
end

"""
    get_cutface_to_iface(m::ImmersedDiscreteModel, d) -> Vector{Int32}

d-cutface → the d-iface it lies on (the inverse of `get_iface_to_cutfaces`).
"""
function get_cutface_to_iface(m::ImmersedDiscreteModel,d::Integer)
  iface_to_cutfaces = get_iface_to_cutfaces(m,d)
  flatten_partition(iface_to_cutfaces,length(iface_to_cutfaces.data))
end

"""
    select_cutfaces(
      m::ImmersedDiscreteModel, d; 
      cutface_mask=nothing, iface_mask=nothing, bgcell_mask=nothing
    ) -> Vector{Int}

Returns the ids of the d-cutfaces that pass the three masks: 
Is in on the cutfaces themselves, on the d-ifaces they lie
on and on the background cells they lie in (`nothing`: no condition).
"""
function select_cutfaces(
  m::ImmersedDiscreteModel,d::Integer;
  cutface_mask=nothing,iface_mask=nothing,bgcell_mask=nothing
)
  cutface_to_bgcell = get_cutfaces(m,d).cell_to_bgcell
  cutface_to_iface = isnothing(iface_mask) ? nothing : get_cutface_to_iface(m,d)
  pass(mask,i) = isnothing(mask) || mask[i]
  findall(eachindex(cutface_to_bgcell)) do c
    pass(cutface_mask,c) &&
      (isnothing(iface_mask) || iface_mask[cutface_to_iface[c]]) &&
      pass(bgcell_mask,cutface_to_bgcell[c])
  end
end

"""
    get_bgcell_to_cutfaces(m::ImmersedDiscreteModel, d) -> Table

Background cell → the d-cutfaces in it (the inverse of the cutface → bgcell glue).
"""
function get_bgcell_to_cutfaces(m::ImmersedDiscreteModel,d::Integer)
  inverse_table(get_cutfaces(m,d).cell_to_bgcell,num_cells(m.bgmodel))
end

"""
    get_bgcell_to_ifaces(m::ImmersedDiscreteModel, d) -> Table

Background cell → the d-ifaces with a cutface in it (sorted).
"""
function get_bgcell_to_ifaces(m::ImmersedDiscreteModel,d::Integer)
  cutface_to_iface = get_cutface_to_iface(m,d)
  cutface_to_ifaces = Table(cutface_to_iface,Int32(1):Int32(length(cutface_to_iface)+1))
  Gridap.Arrays.merge_entries(cutface_to_ifaces,get_bgcell_to_cutfaces(m,d);post=x -> sort!(collect(x)))
end

"""
    SubFaceTriangulation(m::ImmersedDiscreteModel, d)

The d-cutfaces as a triangulation glued to the background model.
"""
SubFaceTriangulation(m::ImmersedDiscreteModel,d::Integer) = SubFaceTriangulation(get_cutfaces(m,d),m.bgmodel)

## As an EmbeddedDiscretization

# The CSG leaf GridapEmbedded addresses the level-set data of an `EmbeddedDiscretization` by
# (as STLCutters' `STLGeometry`): its data is the fitted model.
struct FittedGeometry{T<:DiscreteModel} <: CSG.Geometry
  tree::Leaf{Tuple{T,String,Nothing}}
end
FittedGeometry(fit::DiscreteModel;name="fitted") = FittedGeometry(Leaf((fit,name,nothing)))
CSG.get_tree(geo::FittedGeometry) = geo.tree

"""
    EmbeddedDiscretization(m::ImmersedDiscreteModel)

The cut of the background model of `m` by Γ = ∂Ω, where Ω is the immersed model.
"""
function EmbeddedDiscretization(m::ImmersedDiscreteModel{D}) where D
  fit, bgmodel = get_imodel(m), get_bgmodel(m)

  # Subfacets: the cutfaces of the boundary facets of fit, with their outward normals
  Γ = BoundaryTriangulation(fit)
  bface_to_iface = get_glue(Γ,Val(D-1)).tface_to_mface
  bface_to_normal = collect(lazy_map(evaluate,get_facet_normal(Γ),Fill(zero(Point{D-1,Float64}),num_cells(Γ))))
  iface_to_bface = zeros(Int32,num_faces(fit,D-1))
  iface_to_bface[bface_to_iface] .= 1:length(bface_to_iface)
  is_boundary = Gridap.Geometry.get_isboundary_face(get_grid_topology(fit),D-1)
  facet_ids = select_cutfaces(m,D-1;iface_mask=is_boundary)
  facets = ImmersedFaceData(get_cutfaces(m,D-1),facet_ids)
  normals = bface_to_normal[iface_to_bface[get_cutface_to_iface(m,D-1)[facet_ids]]]
  subfacets = SubFacetData(facets,normals)

  # IN/OUT/CUT: CUT if crossed by Γ, IN if it holds cut cells, OUT otherwise
  is_cut = falses(num_cells(bgmodel))
  is_cut[facets.cell_to_bgcell] .= true
  has_cells = falses(num_cells(bgmodel))
  has_cells[get_cutfaces(m,D).cell_to_bgcell] .= true
  bgcell_to_inoutcut = map((c,h) -> c ? Int8(CUT) : (h ? Int8(IN) : Int8(OUT)),is_cut,has_cells)

  # Subcells: the cut cells in the CUT cells
  cells = ImmersedFaceData(get_cutfaces(m,D),select_cutfaces(m,D;bgcell_mask=is_cut))
  subcells = SubCellData(cells)

  return EmbeddedDiscretization(
    bgmodel,[bgcell_to_inoutcut],
    subcells,[fill(Int8(IN),length(cells.cell_to_bgcell))],
    subfacets,[fill(Int8(INTERFACE),length(facets.cell_to_bgcell))],
    Dict(objectid(fit) => 1),FittedGeometry(fit)
  )
end
