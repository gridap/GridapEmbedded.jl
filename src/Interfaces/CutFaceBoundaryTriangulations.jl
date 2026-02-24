
# Ghost triangulations

function generate_ghost_trian(
  trian::CompositeTriangulation, bgmodel
)
  Dc = num_cell_dims(bgmodel)
  cell_glue = get_glue(trian,Val(Dc))
  return generate_ghost_trian(trian,bgmodel,cell_glue)
end

function generate_ghost_trian(
  trian::CompositeTriangulation, bgmodel, cell_glue::SkeletonPair{<:FaceToFaceGlue}
)
  Dc = num_cell_dims(bgmodel)
  topo = get_grid_topology(bgmodel)
  face_to_cell = get_faces(topo,Dc-1,Dc)
  cell_to_face = get_faces(topo,Dc,Dc-1)

  n_bgfaces = num_faces(bgmodel,Dc-1)
  n_faces = num_cells(trian)
  ghost_faces = zeros(Int32,n_faces)
  p_lcell = ones(Int8,n_bgfaces)
  m_lcell = ones(Int8,n_bgfaces)
  for (i,(p_cell, m_cell)) in enumerate(zip(cell_glue.plus.tface_to_mface,cell_glue.minus.tface_to_mface))
    inter = intersect(cell_to_face[p_cell],cell_to_face[m_cell])
    @assert length(inter) == 1
    face = first(inter)
    ghost_faces[i] = face

    nbors = face_to_cell[ghost_faces[i]]
    p_lcell[face] = findfirst(==(p_cell),nbors)
    m_lcell[face] = findfirst(==(m_cell),nbors)
  end

  plus = BoundaryTriangulation(bgmodel,ghost_faces,p_lcell)
  minus = BoundaryTriangulation(bgmodel,ghost_faces,m_lcell)
  return SkeletonTriangulation(plus,minus)
end

function generate_ghost_trian(
  trian::CompositeTriangulation, bgmodel, cell_glue::FaceToFaceGlue
)
  Dc = num_cell_dims(bgmodel)
  topo = get_grid_topology(bgmodel)
  face_to_cell = get_faces(topo,Dc-1,Dc)
  cell_to_face = get_faces(topo,Dc,Dc-1)
  is_boundary(f) = isone(length(view(face_to_cell,f)))

  n_faces = num_cells(trian)
  ghost_faces = zeros(Int32,n_faces)
  for (i,cell) in enumerate(cell_glue.tface_to_mface)
    faces = filter(is_boundary,view(cell_to_face,cell))
    @assert length(faces) == 1 # TODO: This will break if we are in a corner
    face = first(faces)
    ghost_faces[i] = face
  end

  # NOTE: lcell is always 1 for boundary facets
  return BoundaryTriangulation(bgmodel,ghost_faces)
end

"""
    get_ghost_mask(
      face_trian::SubFacetTriangulation{Df,Dc},
      face_model = get_active_model(face_trian)
    ) where {Df,Dc}

Returns a mask for ghost faces. We define ghost faces as the interfaces between two
different cut facets that are located in different background cells.

The second condition is important: In 3D, some cuts subcells may not be simplices.
In this case, we simplexify the subcell. This creates extra cut interfaces that are
interior to a background cell. These are not considered ghost faces.

- In 2D: Dc = 2, Df = 1 -> Ghost faces have dimension 0 (i.e interface points)
- In 3D: Dc = 3, Df = 2 -> Ghost faces have dimension 1 (i.e interface edges)
"""
function get_ghost_mask(
  face_trian::SubFacetTriangulation{Df,Dc},
  face_model = get_active_model(face_trian)
) where {Df,Dc}
  topo = get_grid_topology(face_model)
  face_to_facets = get_faces(topo,Df-1,Df)

  subfacets = face_trian.subfacets
  facet_to_bgcell = subfacets.facet_to_bgcell

  n_faces = num_faces(topo,Df-1)
  face_is_ghost = zeros(Bool,n_faces)
  for face in 1:n_faces
    facets = view(face_to_facets,face)
    is_boundary = isone(length(facets))
    if !is_boundary
      @assert length(facets) == 2
      bgcells = view(facet_to_bgcell,facets)
      is_ghost = (bgcells[1] != bgcells[2])
      face_is_ghost[face] = is_ghost
    end
  end

  return face_is_ghost
end

"""
    struct CutFaceBoundaryTriangulation{Di,Df,Dp} <: Triangulation{Di,Dp}

Triangulation containing the interfaces between subfacets. We always have dimensions

  - Dc :: Dimension of the background mesh
  - Df = Dc-1 :: Dimension of the cut subfacets
  - Di = Dc-2 :: Dimension of the subfacet interfaces

# Properties

- `face_trian` :: Original SubFacetTriangulation, built on top of the background mesh.
- `face_model` :: Subfacet model. Active model for `face_trian`.
- `face_boundary` :: Triangulation of the interfaces between subfacets. It is glued to the `face_model`.
- `cell_boundary` :: Conceptually the same as `face_boundary`, but it is glued to the
  background mesh cells. Created as a CompositeTriangulation between `face_trian` and `face_boundary`.
- `ghost_boundary` :: Triangulation of the background facets that contain each interface.

The "real" triangulation is `cell_boundary`, but we require the other triangulations to
perform complex changes of domain. Most of the `Triangulation` API is delegated to `cell_boundary`.

## Constructors

    Boundary(face_trian::SubFacetTriangulation)
    Skeleton(face_trian::SubFacetTriangulation)

"""
struct CutFaceBoundaryTriangulation{Di,Df,Dp} <: Triangulation{Di,Dp}
  face_model     :: UnstructuredDiscreteModel{Df,Dp}
  face_trian     :: SubFacetTriangulation{Df,Dp}  # Cut Facet -> BG Cell
  cell_boundary  :: CompositeTriangulation{Di,Dp} # Interface -> BG Cell
  face_boundary  :: BoundaryTriangulation{Di,Dp}  # Interface -> Cut Facet
  ghost_boundary :: BoundaryTriangulation{Df,Dp}  # Ghost Facet -> BG Cell
  interface_sign :: AbstractArray{<:Number}
end

function BoundaryTriangulation(face_trian::SubFacetTriangulation)
  bgmodel = get_background_model(face_trian)
  face_model = get_active_model(face_trian)

  face_boundary  = BoundaryTriangulation(face_model)
  cell_boundary  = CompositeTriangulation(face_trian,face_boundary)
  ghost_boundary = generate_ghost_trian(cell_boundary,bgmodel)
  interface_sign = get_interface_sign(cell_boundary,face_trian,ghost_boundary)

  return CutFaceBoundaryTriangulation(
    face_model,face_trian,cell_boundary,face_boundary,ghost_boundary,interface_sign
  )
end

function get_background_model(t::CutFaceBoundaryTriangulation)
  get_background_model(t.cell_boundary)
end

function get_active_model(t::CutFaceBoundaryTriangulation)
  get_active_model(t.cell_boundary)
end

function get_grid(t::CutFaceBoundaryTriangulation)
  get_grid(t.cell_boundary)
end

# Domain changes

function get_glue(ttrian::CutFaceBoundaryTriangulation{Di,Df,Dp},::Val{D}) where {D,Di,Df,Dp}
  get_glue(ttrian.cell_boundary,Val(D))
end

function is_change_possible(
  strian::SubFacetTriangulation,ttrian::CutFaceBoundaryTriangulation
)
  return strian === ttrian.face_trian
end

function CellData.change_domain(
  a::CellField,ttrian::CutFaceBoundaryTriangulation,tdomain::DomainStyle
)
  strian = get_triangulation(a)
  if strian === ttrian
    # 1) CellField defined on the skeleton
    return change_domain(a,DomainStyle(a),tdomain)
  end

  if is_change_possible(strian,ttrian.cell_boundary)
    # 2) CellField defined on the bgmodel
    b = change_domain(a,ttrian.cell_boundary,tdomain)
  elseif strian === ttrian.face_trian
    # 3) CellField defined on the cut facets
    itrian = Triangulation(ttrian.face_model)
    _a = CellData.similar_cell_field(a,CellData.get_data(a),itrian,DomainStyle(a))
    b = change_domain(_a,ttrian.face_boundary,tdomain)
  else
    @notimplemented
  end
  return CellData.similar_cell_field(b,CellData.get_data(b),ttrian,DomainStyle(b))
end

function CellData.change_domain(
  f::CellData.OperationCellField,ttrian::CutFaceBoundaryTriangulation,tdomain::DomainStyle
)
  args = map(i->change_domain(i,ttrian,tdomain),f.args)
  CellData.OperationCellField(f.op,args...)
end

# Normal vector to the cut facets , n_∂Ω
function get_subfacet_normal_vector(trian::CutFaceBoundaryTriangulation)
  n_∂Ω = get_subfacet_facet_normal(trian.cell_boundary,trian.face_trian)
  return GenericCellField(n_∂Ω,trian,ReferenceDomain())
end

# Normal vector to the ghost facets, n_k
function get_ghost_normal_vector(trian::CutFaceBoundaryTriangulation)
  n = get_ghost_facet_normal(trian.cell_boundary,trian.ghost_boundary)
  return GenericCellField(n,trian,ReferenceDomain())
end

# Orientation of the interface
function get_interface_sign(trian::CutFaceBoundaryTriangulation)
  data = lazy_map(constant_field,trian.interface_sign)
  return GenericCellField(data,trian,ReferenceDomain())
end

# TODO: This is only valid when dealing with linear meshes (where normals are constant over facets).
# If we wanted to use higher-order meshes, we would need to generate the geometric map
# going from the facets to the interfaces.
# However, having a high-order background mesh seems quite silly.
function get_ghost_facet_normal(
  itrian::CompositeTriangulation{Di,Dc}, # Interface -> BG Cell
  gtrian::BoundaryTriangulation{Df,Dc}   # Ghost Facet -> BG Cell
) where {Di,Df,Dc}
  n_g = get_facet_normal(gtrian)
  n_i = lazy_map(evaluate,n_g,Fill(zero(VectorValue{Df,Float64}),num_cells(itrian)))
  return lazy_map(constant_field,n_i)
end

# This one would be fine for higher-order meshes.
function get_subfacet_facet_normal(
  itrian::CompositeTriangulation{Di,Dc}, # Interface -> BG Cell
  ftrian::SubFacetTriangulation{Df,Dc},  # Cut Facet -> BG Cell
) where {Di,Df,Dc}
  glue = get_glue(itrian.dtrian,Val(Df))
  i_to_f_ids = glue.tface_to_mface
  i_to_f_map = glue.tface_to_mface_map
  n_f = lazy_map(Reindex(get_facet_normal(ftrian)),i_to_f_ids)
  n_i = lazy_map(Broadcasting(∘),n_f,i_to_f_map)
  return n_i
end

# There is still something sweaty about this...
# Why do we apply the sign change but at the same time call `get_edge_tangents` in the
# creation of conormal vectors in 3D? It's like we are cancelling the sign change...
# There is more to think about here.
function get_interface_sign(
  itrian::CompositeTriangulation{Di,Dc}, # Interface -> BG Cell
  ftrian::SubFacetTriangulation{Df,Dc},  # Cut Facet -> BG Cell
  gtrian::BoundaryTriangulation{Df,Dc},  # Ghost Facet -> BG Cell
) where {Di,Df,Dc}
  function signdot(a,b)
    s = sign(dot(a,b))
    return ifelse(iszero(s),1,s)
  end
  n_∂Ω = get_subfacet_facet_normal(itrian,ftrian)
  n_k = get_ghost_facet_normal(itrian,gtrian)
  if Di == 0
    cross2D(n) = VectorValue(-n[2],n[1])
    n_S = lazy_map(Operation(cross2D),n_k)
  else
    t_S = get_edge_tangents(itrian.dtrian)
    n_S = lazy_map(Operation(cross),n_k,t_S)
  end
  sgn = lazy_map(Operation(signdot),n_∂Ω,n_S)
  return collect(lazy_map(evaluate,sgn,Fill(zero(VectorValue{Di,Float64}),num_cells(itrian))))
end

function get_edge_tangents(trian::BoundaryTriangulation{1})
  function t(c)
    @assert length(c) == 2
    t = c[2] - c[1]
    return t/norm(t)
  end
  return lazy_map(constant_field,lazy_map(t,get_cell_coordinates(trian)))
end

function get_edge_tangents(trian::CutFaceBoundaryTriangulation{1})
  data = get_edge_tangents(trian.face_boundary)
  return GenericCellField(data,trian,ReferenceDomain())
end

# Normal vector to the cut interface, n_S
function get_normal_vector(trian::CutFaceBoundaryTriangulation{Di}) where {Di}
  n_k = get_ghost_normal_vector(trian)
  isign = get_interface_sign(trian)
  if Di == 0 # 2D
    cross2D(n) = VectorValue(-n[2],n[1])
    n_S = Operation(cross2D)(n_k) # nS = nk x tS and tS = ±e₃ in 2D
  elseif Di == 1 # 3D
    t_S = get_edge_tangents(trian)
    n_S = Operation(cross)(n_k,t_S) # nk = tS x nS -> nS = nk x tS (eq 6.25)
  else
    @notimplemented
  end
  return n_S * isign
end

get_facet_normal(trian::CutFaceBoundaryTriangulation) = get_data(get_normal_vector(trian))

# Tangent vector to the cut interface, t_S = n_S x n_k
function get_tangent_vector(trian::CutFaceBoundaryTriangulation{Di}) where {Di}
  @notimplementedif Di != 1
  n_S = get_normal_vector(trian)
  n_k = get_ghost_normal_vector(trian)
  return Operation(cross)(n_S,n_k)
end

# Conormal vectors, m_k = t_S x n_∂Ω
function get_conormal_vector(trian::CutFaceBoundaryTriangulation{Di}) where {Di}
  n_∂Ω = get_subfacet_normal_vector(trian)
  isign = get_interface_sign(trian)
  if Di == 0 # 2D
    cross2D(n) = VectorValue(n[2],-n[1])
    m_k = Operation(cross2D)(n_∂Ω)
  elseif Di == 1 # 3D
    t_S = get_edge_tangents(trian)
    m_k = Operation(cross)(t_S,n_∂Ω) # m_k = t_S x n_∂Ω (eq 6.26)
  else
    @notimplemented
  end
  return m_k * isign
end

# CutFaceSkeletonTriangulation & CutFaceBoundaryTriangulationView
const CutFaceBoundaryTriangulationView{Di,Df,Dp} = TriangulationView{Di,Dp,CutFaceBoundaryTriangulation{Di,Df,Dp}}
const CutFaceSkeletonTriangulation{Di,Df,Dp} = SkeletonTriangulation{Di,Dp,<:Union{
  CutFaceBoundaryTriangulation{Di,Df,Dp},
  CutFaceBoundaryTriangulationView{Di,Df,Dp}
  }
}

function SkeletonTriangulation(face_trian::SubFacetTriangulation)
  bgmodel = get_background_model(face_trian)
  face_model = get_active_model(face_trian)

  ghost_mask = get_ghost_mask(face_trian,face_model)
  face_skeleton = SkeletonTriangulation(face_model,ghost_mask)
  cell_skeleton = CompositeTriangulation(face_trian,face_skeleton)
  ghost_skeleton = generate_ghost_trian(cell_skeleton,bgmodel)

  ctrian_plus = CompositeTriangulation(face_trian,face_skeleton.plus)
  ctrian_minus = CompositeTriangulation(face_trian,face_skeleton.minus)
  isign_plus = get_interface_sign(ctrian_plus,face_trian,ghost_skeleton.plus)
  isign_minus = get_interface_sign(ctrian_plus,face_trian,ghost_skeleton.minus)

  plus = CutFaceBoundaryTriangulation(
    face_model,face_trian,ctrian_plus,
    face_skeleton.plus,ghost_skeleton.plus,isign_plus
  )
  minus = CutFaceBoundaryTriangulation(
    face_model,face_trian,ctrian_minus,
    face_skeleton.minus,ghost_skeleton.minus,isign_minus
  )
  return SkeletonTriangulation(plus,minus)
end

for func in (:get_subfacet_normal_vector,:get_ghost_normal_vector,:get_conormal_vector)
  @eval begin
    function $func(trian::CutFaceSkeletonTriangulation)
      plus  = GenericCellField(CellData.get_data($func(trian.plus)),trian,ReferenceDomain())
      minus = GenericCellField(CellData.get_data($func(trian.minus)),trian,ReferenceDomain())
      return SkeletonPair(plus,minus)
    end
  end
end

for func in (:get_normal_vector,:get_tangent_vector)
  @eval begin
    function $func(trian::CutFaceSkeletonTriangulation)
      return GenericCellField(CellData.get_data($func(trian.plus)),trian,ReferenceDomain())
    end
  end
end

for func in (:get_tangent_vector,:get_subfacet_normal_vector,:get_ghost_normal_vector,:get_conormal_vector)
  @eval begin
    function $func(trian::CutFaceBoundaryTriangulationView)
      data = CellData.get_data($func(trian.parent))
      restricted_data = restrict(data,trian.cell_to_parent_cell)
      return GenericCellField(restricted_data,trian,ReferenceDomain())
    end
  end
end

############################################################################################
# This will go to Gridap
#
# function Arrays.evaluate!(cache,k::Operation,a::SkeletonPair{<:CellField})
#   plus = k(a.plus)
#   minus = k(a.minus)
#   SkeletonPair(plus,minus)
# end

# function Arrays.evaluate!(cache,k::Operation,a::SkeletonPair{<:CellField},b::SkeletonPair{<:CellField})
#   plus = k(a.plus,b.plus)
#   minus = k(a.minus,b.minus)
#   SkeletonPair(plus,minus)
# end

# import Gridap.TensorValues: inner, outer
# import LinearAlgebra: dot
# import Base: abs, *, +, -, /

# for op in (:/,)
#   @eval begin
#     ($op)(a::CellField,b::SkeletonPair{<:CellField}) = Operation($op)(a,b)
#     ($op)(a::SkeletonPair{<:CellField},b::CellField) = Operation($op)(a,b)
#   end
# end

# for op in (:outer,:*,:dot,:/)
#   @eval begin
#     ($op)(a::SkeletonPair{<:CellField},b::SkeletonPair{<:CellField}) = Operation($op)(a,b)
#   end
# end

# function CellData.change_domain(a::SkeletonPair, ::ReferenceDomain, ::PhysicalDomain)
#   plus = change_domain(a.plus,ReferenceDomain(),PhysicalDomain())
#   minus = change_domain(a.minus,ReferenceDomain(),PhysicalDomain())
#   return SkeletonPair(plus,minus)
# end
