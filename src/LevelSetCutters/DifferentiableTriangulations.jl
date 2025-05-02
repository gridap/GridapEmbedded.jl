

"""
    mutable struct DifferentiableTriangulation{Dc,Dp} <: Triangulation{Dc,Dp}

A DifferentiableTriangulation is a wrapper around an embedded triangulation
(i.e SubCellTriangulation or SubFacetTriangulation) implementing all the necessary
methods to compute derivatives w.r.t. deformations of the embedded mesh.

To do so, it propagates dual numbers into the geometric maps mapping cut subcells/subfacets
to the background mesh.

## Constructor: 

    DifferentiableTriangulation(trian::Triangulation,fe_space::FESpace)

where `trian` must be an embedded triangulation and `fe_space` is the `FESpace` where 
the level-set function lives.

"""
mutable struct DifferentiableTriangulation{Dc,Dp,A,B} <: Triangulation{Dc,Dp}
  trian :: A
  fe_space :: B
  cell_values
  caches
  function DifferentiableTriangulation(
    trian :: Triangulation{Dc,Dp},
    fe_space :: FESpace,
    cell_values,caches
  ) where {Dc,Dp}
    A = typeof(trian)
    B = typeof(fe_space)
    new{Dc,Dp,A,B}(trian,fe_space,cell_values,caches)
  end
end

# Constructors

DifferentiableTriangulation(trian::Triangulation,fe_space) = trian

function DifferentiableTriangulation(
  trian::Union{<:SubCellTriangulation,<:SubFacetTriangulation},
  fe_space::FESpace
)
  caches = precompute_autodiff_caches(trian)
  return DifferentiableTriangulation(trian,fe_space,nothing,caches)
end

# Update cell values

(t::DifferentiableTriangulation)(φh) = update_trian!(t,get_fe_space(φh),φh)

update_trian!(trian::Triangulation,U,φh) = trian

function update_trian!(trian::DifferentiableTriangulation,space::FESpace,φh)
  (trian.fe_space !== space) && return trian
  trian.cell_values = extract_dualized_cell_values(trian.trian,φh)
  return trian
end

function update_trian!(trian::DifferentiableTriangulation,::FESpace,::Nothing)
  trian.cell_values = nothing
  return trian
end

# Autodiff

function FESpaces._change_argument(
  op,f,trian::DifferentiableTriangulation,uh
)
  U = get_fe_space(uh)
  function g(cell_u)
    cf = CellField(U,cell_u)
    update_trian!(trian,U,cf)
    cell_grad = f(cf)
    update_trian!(trian,U,nothing) # TODO: experimental
    get_contribution(cell_grad,trian)
  end
  g
end

function FESpaces._compute_cell_ids(uh,ttrian::DifferentiableTriangulation)
  FESpaces._compute_cell_ids(uh,ttrian.trian)
end

function Geometry.get_background_model(t::DifferentiableTriangulation)
  get_background_model(t.trian)
end

function Geometry.get_grid(t::DifferentiableTriangulation)
  get_grid(t.trian)
end

function Geometry.get_cell_reffe(t::DifferentiableTriangulation)
  get_cell_reffe(t.trian)
end

# TODO: Do we ever need to dualize the cell points?
# I think its not necessary, since all the dual numbers are propagated through the cellmaps...
# Also: The current version dualizes only the phys points...
# If we want to indeed dualize this, we should probably also dualize the ref points
# in the case where ttrian.trian is a SubCellTriangulation (but not in the case of a SubFacetTriangulation)
# Anyway, I don't think this matters for now...
function CellData.get_cell_points(ttrian::DifferentiableTriangulation)
  pts = get_cell_points(ttrian.trian)
  cell_ref_point = pts.cell_ref_point
  if isnothing(ttrian.cell_values) || isempty(ttrian.cell_values)
    cell_phys_point = pts.cell_phys_point
  else
    c = ttrian.caches
    cell_phys_point = lazy_map(
      DualizeCoordsMap(),c.face_to_coords,c.face_to_bgcoords,
      ttrian.cell_values,c.face_to_edges,c.face_to_edge_lists
    )
  end
  return CellPoint(cell_ref_point, cell_phys_point, ttrian, DomainStyle(pts))
end

function Geometry.get_cell_map(ttrian::DifferentiableTriangulation)
  if isnothing(ttrian.cell_values) || isempty(ttrian.cell_values)
    return get_cell_map(ttrian.trian)
  end
  c = ttrian.caches
  cell_values = ttrian.cell_values
  cell_to_coords = lazy_map(
    DualizeCoordsMap(),c.face_to_coords,c.face_to_bgcoords,
    cell_values,c.face_to_edges,c.face_to_edge_lists
  )
  cell_reffe = get_cell_reffe(ttrian)
  cell_map = compute_cell_maps(cell_to_coords,cell_reffe)
  return cell_map
end

function face_normal(face_coords::Vector{<:Point},orientation::Int8)
  n = face_normal(face_coords)
  return n*orientation
end
function face_normal(face_coords::Vector{<:Point{2}})
  p1, p2 = face_coords[1:2]
  LevelSetCutters._normal_vector(p2-p1)
end
function face_normal(face_coords::Vector{<:Point{3}})
  p1, p2, p3 = face_coords[1:3]
  LevelSetCutters._normal_vector(p2-p1,p3-p1)
end

function Geometry.get_facet_normal(
  ttrian::DifferentiableTriangulation{Dc,Dp,<:SubFacetTriangulation}
) where {Dc,Dp}
  if isnothing(ttrian.cell_values) || isempty(ttrian.cell_values)
    return get_facet_normal(ttrian.trian)
  end
  c = ttrian.caches
  cell_values = ttrian.cell_values
  cell_to_coords = lazy_map(
    DualizeCoordsMap(),c.face_to_coords,c.face_to_bgcoords,
    cell_values,c.face_to_edges,c.face_to_edge_lists
  )
  facet_normals = lazy_map(face_normal,cell_to_coords,c.orientations)
  return lazy_map(constant_field,facet_normals)
end

function Geometry.get_glue(ttrian::DifferentiableTriangulation,val::Val{D}) where {D}
  glue = get_glue(ttrian.trian,val)
  if isnothing(glue) || isnothing(ttrian.cell_values) || isempty(ttrian.cell_values)
    return glue
  end

  # New reference maps
  c = ttrian.caches
  cell_values = ttrian.cell_values
  cell_to_rcoords = lazy_map(
    DualizeCoordsMap(),c.face_to_rcoords,c.face_to_bgrcoords,
    cell_values,c.face_to_edges,c.face_to_edge_lists
  )
  cell_reffe = get_cell_reffe(ttrian)
  ref_cell_map = compute_cell_maps(cell_to_rcoords,cell_reffe)

  return FaceToFaceGlue(
    glue.tface_to_mface,
    ref_cell_map,
    glue.mface_to_tface,
  )
end

function Geometry.is_change_possible(
  strian::A,ttrian::DifferentiableTriangulation{Dc,Dp,A}
) where {Dc,Dp,A <: Union{SubCellTriangulation,SubFacetTriangulation}}
  return strian === ttrian.trian
end

function Geometry.best_target(
  strian::A,ttrian::DifferentiableTriangulation{Dc,Dp,A}
) where {Dc,Dp,A <: Union{SubCellTriangulation,SubFacetTriangulation}}
  return ttrian
end

for tdomain in (:ReferenceDomain,:PhysicalDomain)
  for sdomain in (:ReferenceDomain,:PhysicalDomain)
    @eval begin
      function CellData.change_domain(
        a::CellField,strian::A,::$sdomain,ttrian::DifferentiableTriangulation{Dc,Dp,A},::$tdomain
      ) where {Dc,Dp,A <: Union{SubCellTriangulation,SubFacetTriangulation}}
        @assert is_change_possible(strian,ttrian)
        b = change_domain(a,$(tdomain)())
        return CellData.similar_cell_field(a,CellData.get_data(b),ttrian,$(tdomain)())
      end
    end
  end
end

function FESpaces.get_cell_fe_data(fun,f,ttrian::DifferentiableTriangulation)
  FESpaces.get_cell_fe_data(fun,f,ttrian.trian)
end

function compute_cell_maps(cell_coords,cell_reffes)
  cell_shapefuns = lazy_map(get_shapefuns,cell_reffes)
  default_cell_map = lazy_map(linear_combination,cell_coords,cell_shapefuns)
  default_cell_grad = lazy_map(∇,default_cell_map)
  cell_poly = lazy_map(get_polytope,cell_reffes)
  cell_q0 = lazy_map(p->zero(first(get_vertex_coordinates(p))),cell_poly)
  origins = lazy_map(evaluate,default_cell_map,cell_q0)
  gradients = lazy_map(evaluate,default_cell_grad,cell_q0)
  cell_map = lazy_map(Fields.affine_map,gradients,origins)
  return cell_map
end

# DualizeCoordsMap

struct DualizeCoordsMap <: Map end

function Arrays.return_cache(
  k::DualizeCoordsMap,
  coords::Vector{<:Point{Dp,Tp}},
  bg_coords::Vector{<:Point{Dp,Tp}},
  values::Vector{Tv},
  edges::Vector{Int8},
  edge_list::Vector{Vector{Int8}}
) where {Dp,Tp,Tv}
  T = Point{Dp,Tv}
  return CachedArray(zeros(T, length(coords)))
end

function Arrays.evaluate!(
  cache,
  k::DualizeCoordsMap,
  coords::Vector{<:Point{Dp,Tp}},
  bg_coords::Vector{<:Point{Dp,Tp}},
  values::Vector{Tv},
  edges::Vector{Int8},
  edge_list::Vector{Vector{Int8}}
) where {Dp,Tp,Tv}
  setsize!(cache,(length(coords),))
  new_coords = cache.array
  for (i,e) in enumerate(edges)
    if e == -1
      new_coords[i] = coords[i]
    else
      n1, n2 = edge_list[e]
      q1, q2 = bg_coords[n1], bg_coords[n2]
      v1, v2 = abs(values[n1]), abs(values[n2])
      λ = v1/(v1+v2)
      new_coords[i] = q1 + λ*(q2-q1)
    end
  end
  return new_coords
end

"""
    precompute_cut_edge_ids(rcoords,bg_rcoords,edge_list)

Given
  - `rcoords`: the node ref coordinates of the cut subcell/subfacet,
  - `bg_rcoords`: the node ref coordinates of the background cell containing it,
  - `edge_list`: the list of nodes defining each edge of the background cell,

this function returns a vector that for each node of the cut subcell/subfacet contains
  - `-1` if the node is also a node of the background cell,
  - the id of the edge containing the node otherwise.
"""
function precompute_cut_edge_ids(
  rcoords::Vector{<:Point{Dp,Tp}},
  bg_rcoords::Vector{<:Point{Dp,Tp}},
  edge_list::Vector{<:Vector{<:Integer}}
) where {Dp,Tp}
  tol = 10*eps(Tp)
  edges = Vector{Int8}(undef,length(rcoords))
  for (i,p) in enumerate(rcoords)
    if any(q -> norm(q-p) < tol, bg_rcoords)
      edges[i] = Int8(-1)
    else
      e = findfirst(edge -> belongs_to_edge(p,edge,bg_rcoords), edge_list)
      edges[i] = Int8(e)
    end
  end
  return edges
end

function get_edge_list(poly::Polytope)
  ltcell_to_lpoints, simplex = simplexify(poly)
  simplex_edges = get_faces(simplex,1,0)
  ltcell_to_edges = map(pts -> map(e -> pts[e], simplex_edges), ltcell_to_lpoints)
  return collect(Vector{Int8},unique(sort,vcat(ltcell_to_edges...)))
end

function belongs_to_edge(
  p::Point{D,T},edge::Vector{<:Integer},bgpts::Vector{Point{D,T}}
) where {D,T}
  tol = 10*eps(T)
  p1, p2 = bgpts[edge]
  return norm(cross(p-p1,p2-p1)) < tol
end

function precompute_autodiff_caches(
  trian::SubCellTriangulation
)
  bgmodel = get_background_model(trian)
  subcells = trian.subcells

  precompute_autodiff_caches(
    bgmodel,
    subcells.cell_to_bgcell,
    subcells.cell_to_points,
    subcells.point_to_rcoords,
    subcells.point_to_coords,
  )
end

function precompute_autodiff_caches(
  trian::SubFacetTriangulation
)
  bgmodel = get_background_model(trian)
  subfacets = trian.subfacets

  caches = precompute_autodiff_caches(
    bgmodel,
    subfacets.facet_to_bgcell,
    subfacets.facet_to_points,
    subfacets.point_to_rcoords,
    subfacets.point_to_coords,
  )

  # Precompute orientations
  orientations = collect(lazy_map(orient,subfacets.facet_to_normal,caches.face_to_coords))

  cache = (; caches..., orientations)
  return cache
end

orient(n,fcoords) = round(Int8,dot(n,face_normal(fcoords)))
Arrays.return_value(::typeof(orient),n,face_coords) = zero(Int8)

function precompute_autodiff_caches(
  bgmodel,
  face_to_bgcell,
  face_to_points,
  point_to_rcoords,
  point_to_coords,
)
  bg_ctypes = get_cell_type(bgmodel)
  bgcell_to_polys = expand_cell_data(get_polytopes(bgmodel),bg_ctypes)
  bgcell_to_coords = get_cell_coordinates(bgmodel)
  bgcell_to_rcoords = lazy_map(get_vertex_coordinates,bgcell_to_polys)

  face_to_bgcoords = lazy_map(Reindex(bgcell_to_coords),face_to_bgcell)
  face_to_bgrcoords = lazy_map(Reindex(bgcell_to_rcoords),face_to_bgcell)
  face_to_rcoords = lazy_map(Broadcasting(Reindex(point_to_rcoords)),face_to_points)
  face_to_coords = lazy_map(Broadcasting(Reindex(point_to_coords)),face_to_points)

  bgcell_to_edge_lists = lazy_map(get_edge_list,bgcell_to_polys)
  face_to_edge_lists = lazy_map(Reindex(bgcell_to_edge_lists),face_to_bgcell)
  face_to_edges = collect(lazy_map(precompute_cut_edge_ids,face_to_rcoords,face_to_bgrcoords,face_to_edge_lists))

  cache = (;
    face_to_rcoords,
    face_to_coords,
    face_to_bgrcoords,
    face_to_bgcoords,
    face_to_edges,
    face_to_edge_lists
  )
  return cache
end

function extract_dualized_cell_values(
  trian::SubCellTriangulation,
  φh::CellField,
)
  @assert isa(DomainStyle(φh),ReferenceDomain)
  bgmodel = get_background_model(trian)
  bgcell_to_values = extract_dualized_cell_values(bgmodel,φh)

  subcells = trian.subcells
  cell_to_bgcell = subcells.cell_to_bgcell
  cell_to_values = lazy_map(Reindex(bgcell_to_values),cell_to_bgcell)
  return cell_to_values
end

function extract_dualized_cell_values(
  trian::SubFacetTriangulation,
  φh::CellField,
)
  @assert isa(DomainStyle(φh),ReferenceDomain)
  bgmodel = get_background_model(trian)
  bgcell_to_values = extract_dualized_cell_values(bgmodel,φh)

  subfacets = trian.subfacets
  facet_to_bgcell = subfacets.facet_to_bgcell
  facet_to_values = lazy_map(Reindex(bgcell_to_values),facet_to_bgcell)
  return facet_to_values
end

function extract_dualized_cell_values(
  bgmodel::DiscreteModel,
  φh::CellField,
)
  @assert isa(DomainStyle(φh),ReferenceDomain)
  bg_ctypes = get_cell_type(bgmodel)
  bgcell_to_polys = expand_cell_data(get_polytopes(bgmodel),bg_ctypes)
  bgcell_to_rcoords = lazy_map(get_vertex_coordinates,bgcell_to_polys)
  bgcell_to_fields = CellData.get_data(φh)
  bgcell_to_values = lazy_map(evaluate,bgcell_to_fields,bgcell_to_rcoords)
  return bgcell_to_values
end

# TriangulationView
# This is mostly used in distributed, where we remove ghost cells by taking a view 
# of the local triangulations. 

const DifferentiableTriangulationView{Dc,Dp} = Geometry.TriangulationView{Dc,Dp,<:DifferentiableTriangulation}

function DifferentiableTriangulation(
  trian :: Geometry.TriangulationView,
  fe_space :: FESpace
)
  parent = DifferentiableTriangulation(trian.parent,fe_space)
  return Geometry.TriangulationView(parent,trian.cell_to_parent_cell)
end

function update_trian!(trian::Geometry.TriangulationView,U,φh)
  update_trian!(trian.parent,U,φh)
  return trian
end

function FESpaces._change_argument(
  op,f,trian::DifferentiableTriangulationView,uh
)
  U = get_fe_space(uh)
  function g(cell_u)
    cf = CellField(U,cell_u)
    update_trian!(trian,U,cf)
    cell_grad = f(cf)
    update_trian!(trian,U,nothing)
    get_contribution(cell_grad,trian)
  end
  g
end

# AppendedTriangulation
#
# When cutting an embedded domain, we will usually end up with an AppendedTriangulation
# containing
#   a) a regular triangulation with the IN/OUT cells
#   b) a SubCell/SubFacetTriangulation with the CUT cells
# We only need to propagate the dual numbers to the CUT cells, which is what the
# following implementation does:

const DifferentiableAppendedTriangulation{Dc,Dp,A} = 
  AppendedTriangulation{Dc,Dp,<:Union{<:DifferentiableTriangulation,<:DifferentiableTriangulationView{Dc,Dp}}}

function DifferentiableTriangulation(
  trian::AppendedTriangulation, fe_space::FESpace
)
  a = DifferentiableTriangulation(trian.a,fe_space)
  b = DifferentiableTriangulation(trian.b,fe_space)
  return AppendedTriangulation(a,b)
end

function update_trian!(trian::DifferentiableAppendedTriangulation,U,φh)
  update_trian!(trian.a,U,φh)
  update_trian!(trian.b,U,φh)
  return trian
end

function FESpaces._change_argument(
  op,f,trian::DifferentiableAppendedTriangulation,uh
)
  U = get_fe_space(uh)
  function g(cell_u)
    cf = CellField(U,cell_u)
    update_trian!(trian,U,cf)
    cell_grad = f(cf)
    update_trian!(trian,U,nothing)
    get_contribution(cell_grad,trian)
  end
  g
end

# TODO: Move to Gridap
function FESpaces._compute_cell_ids(uh,ttrian::AppendedTriangulation)
  ids_a = FESpaces._compute_cell_ids(uh,ttrian.a)
  ids_b = FESpaces._compute_cell_ids(uh,ttrian.b)
  lazy_append(ids_a,ids_b)
end
