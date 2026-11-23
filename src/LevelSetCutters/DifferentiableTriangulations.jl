

"""
    mutable struct DifferentiableTriangulation{Dc,Dp} <: Triangulation{Dc,Dp}

A DifferentiableTriangulation is a wrapper around an embedded triangulation
(i.e a SubFaceTriangulation: SubCellTriangulation or SubFacetTriangulation) implementing all the necessary
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

const EmbeddedTriangulationTypes = Union{SubFaceTriangulation,SubFacetBoundaryTriangulation}

# Constructors

DifferentiableTriangulation(trian::Triangulation,fe_space) = trian

function DifferentiableTriangulation(
  trian::SubFaceTriangulation,
  fe_space::FESpace
)
  caches = precompute_autodiff_caches(trian)
  return DifferentiableTriangulation(trian,fe_space,nothing,caches)
end

function DifferentiableTriangulation(
  trian::SubFacetBoundaryTriangulation,
  fe_space::FESpace
)
  caches = precompute_autodiff_caches(trian,fe_space)
  return DifferentiableTriangulation(trian,fe_space,nothing,caches)
end

# Update cell values

(t::DifferentiableTriangulation)(φh) = update_trian!(t,get_fe_space(φh),φh)

update_trian!(trian::Triangulation,U,φh) = trian

function update_trian!(trian::DifferentiableTriangulation,space::FESpace,φh)
  (trian.fe_space !== space) && return trian
  trian.cell_values = extract_dualized_cell_values(trian.trian,φh,trian.caches)
  return trian
end

function update_trian!(trian::DifferentiableTriangulation,::FESpace,::Nothing)
  trian.cell_values = nothing
  return trian
end

const MultiFieldSpaceTypes = Union{<:MultiFieldFESpace,<:DistributedMultiFieldFESpace}

function update_trian!(trian::DifferentiableTriangulation,space::MultiFieldSpaceTypes,φh)
  map((Ui,φi)->update_trian!(trian,Ui,φi),space,φh)
  return trian
end

function update_trian!(trian::DifferentiableTriangulation,::MultiFieldSpaceTypes,::Nothing)
  trian.cell_values = nothing
  return trian
end

# Autodiff

_is_levelset_space(f::FESpace,ttrian::DifferentiableTriangulation) = (f === ttrian.fe_space)
_is_levelset_space(f::FEFunction,ttrian::DifferentiableTriangulation) = (get_fe_space(f) === ttrian.fe_space)
_is_levelset_space(f,ttrian::DifferentiableTriangulation) = false

function FESpaces._compute_cell_ids(uh,ttrian::DifferentiableTriangulation)
  FESpaces._compute_cell_ids(uh,ttrian.trian)
end

function FESpaces._compute_cell_ids(
  uh,ttrian::DifferentiableTriangulation{Dc,Dp,<:SubFacetBoundaryTriangulation}
) where {Dc,Dp}
  _is_levelset_space(uh,ttrian) ? ttrian.caches.face_to_φscell : FESpaces._compute_cell_ids(uh,ttrian.trian)
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

Geometry.get_facet_normal(ttrian::DifferentiableTriangulation) = get_facet_normal(ttrian.trian)

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

function _dualized_ref_cell_map(ttrian::DifferentiableTriangulation)
  c = ttrian.caches
  cell_to_rcoords = lazy_map(
    DualizeCoordsMap(),c.face_to_rcoords,c.face_to_bgrcoords,
    ttrian.cell_values,c.face_to_edges,c.face_to_edge_lists
  )
  compute_cell_maps(cell_to_rcoords,get_cell_reffe(ttrian))
end

function Geometry.get_glue(ttrian::DifferentiableTriangulation,val::Val{D}) where {D}
  glue = get_glue(ttrian.trian,val)
  (isnothing(glue) || isnothing(ttrian.cell_values) || isempty(ttrian.cell_values)) && return glue
  FaceToFaceGlue(glue.tface_to_mface,_dualized_ref_cell_map(ttrian),glue.mface_to_tface)
end

function Geometry.get_glue(
  ttrian::DifferentiableTriangulation{Dc,Dp,<:SubFacetBoundaryTriangulation},val::Val{D}
) where {Dc,Dp,D}
  glue = get_glue(ttrian.trian,val)
  (isnothing(glue) || isnothing(ttrian.cell_values) || isempty(ttrian.cell_values)) && return glue
  subfacet_to_facet_map = _dualized_ref_cell_map(ttrian)
  if D == Dc
    FaceToFaceGlue(glue.tface_to_mface,subfacet_to_facet_map,glue.mface_to_tface)
  elseif D - 1 == Dc
    fglue = get_glue(ttrian.trian.facets,val)
    face_ref_map = lazy_map(Reindex(fglue.tface_to_mface_map),ttrian.trian.subfacet_to_facet)
    cell_ref_map = lazy_map(∘,face_ref_map,subfacet_to_facet_map)
    FaceToFaceGlue(glue.tface_to_mface,cell_ref_map,glue.mface_to_tface)
  else
    nothing
  end
end

function Geometry.is_change_possible(
  strian::A,ttrian::DifferentiableTriangulation{Dc,Dp,A}
) where {Dc,Dp,A <: EmbeddedTriangulationTypes}
  return strian === ttrian.trian
end

function Geometry.best_target(
  strian::A,ttrian::DifferentiableTriangulation{Dc,Dp,A}
) where {Dc,Dp,A <: EmbeddedTriangulationTypes}
  return ttrian
end

for tdomain in (:ReferenceDomain,:PhysicalDomain)
  for sdomain in (:ReferenceDomain,:PhysicalDomain)
    @eval begin
      function CellData.change_domain(
        a::CellField,strian::A,::$sdomain,ttrian::DifferentiableTriangulation{Dc,Dp,A},::$tdomain
      ) where {Dc,Dp,A <: EmbeddedTriangulationTypes}
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

function FESpaces.get_cell_fe_data(
  fun,f,ttrian::DifferentiableTriangulation{Dc,Dp,<:SubFacetBoundaryTriangulation}
) where {Dc,Dp}
  _is_levelset_space(f,ttrian) ?
    lazy_map(Reindex(fun(f)),ttrian.caches.face_to_φscell) :
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
  p::Point{1,T},edge::Vector{<:Integer},bgpts::Vector{Point{1,T}}
) where {T}
  tol = 10*eps(T)
  a, b = minmax(bgpts[edge[1]][1],bgpts[edge[2]][1])
  return a - tol <= p[1] <= b + tol
end

function belongs_to_edge(
  p::Point{D,T},edge::Vector{<:Integer},bgpts::Vector{Point{D,T}}
) where {D,T}
  tol = 10*eps(T)
  p1, p2 = bgpts[edge]
  return norm(cross(p-p1,p2-p1)) < tol
end

function precompute_autodiff_caches(
  trian::SubFaceTriangulation
)
  bgmodel = get_background_model(trian)
  subfaces = trian.subfaces

  caches = precompute_autodiff_caches(
    bgmodel,
    subfaces.cell_to_bgcell,
    subfaces.cell_to_points,
    subfaces.point_to_rcoords,
    subfaces.point_to_coords,
  )
  cell_to_normal = subfaces.cell_to_normal
  isnothing(cell_to_normal) && return caches

  # Precompute orientations (stored normals vs normals from the vertex order)
  orientations = collect(lazy_map(orient,cell_to_normal,caches.face_to_coords))

  cache = (; caches..., orientations)
  return cache
end

function precompute_autodiff_caches(
  trian::SubFacetBoundaryTriangulation{Dc,Dp},
  fe_space::FESpace
) where {Dc,Dp}
  bgmodel = get_background_model(trian)
  subfacets = trian.subfacets
  bgfacet_grid = Grid(ReferenceFE{Dc},bgmodel)
  caches = precompute_autodiff_caches(
    bgfacet_grid,
    subfacets.cell_to_bgcell,
    subfacets.cell_to_points,
    subfacets.point_to_rcoords,
    subfacets.point_to_coords,
  )

  # Find the adjacent cell K_φ ∈ get_triangulation(fe_space) for each cut background facet
  D = num_cell_dims(bgmodel)
  φtrian = get_triangulation(fe_space)
  @check get_background_model(φtrian) === bgmodel
  φglue = get_glue(φtrian,Val(D))
  @check isa(φglue,FaceToFaceGlue) && !isnothing(φglue.mface_to_tface)
  bgcell_to_φcell = φglue.mface_to_tface

  topo = get_grid_topology(bgmodel)
  facet_to_cells = Table(get_faces(topo,Dc,D))
  bgfacets = subfacets.cell_to_bgcell
  n_subfacets = length(bgfacets)
  face_to_φscell = Vector{Int32}(undef,n_subfacets)
  face_to_φlcell = Vector{Int8}(undef,n_subfacets)
  for sf in 1:n_subfacets
    bgfacet = bgfacets[sf]
    cells = view(facet_to_cells,bgfacet)
    eval_cell = trian.cell_ids[sf]
    idx = findfirst(c -> c == eval_cell && bgcell_to_φcell[c] > 0, cells)
    isnothing(idx) && (idx = findfirst(c -> bgcell_to_φcell[c] > 0, cells))
    @check !isnothing(idx) "Cut background facet $bgfacet has no adjacent cell in the level-set FESpace."
    face_to_φscell[sf] = Int32(bgcell_to_φcell[cells[idx]])
    face_to_φlcell[sf] = Int8(idx)
  end

  cell_grid = get_grid(bgmodel)
  φ_face_grid = view(bgfacet_grid,bgfacets)
  φ_glue = Geometry.OverlappingFaceToCellGlue(topo,cell_grid,φ_face_grid,bgfacets,face_to_φlcell)
  face_to_φrcoords = collect(Geometry._compute_face_to_q_vertex_coords(cell_grid,φ_face_grid,φ_glue))

  return (; caches..., face_to_φscell, face_to_φrcoords)
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

extract_dualized_cell_values(trian,φh,caches) = extract_dualized_cell_values(trian,φh)

function extract_dualized_cell_values(
  trian::SubFaceTriangulation,
  φh::CellField,
)
  @assert isa(DomainStyle(φh),ReferenceDomain)
  bgmodel = get_background_model(trian)
  bgcell_to_values = extract_dualized_cell_values(bgmodel,φh)

  cell_to_bgcell = trian.subfaces.cell_to_bgcell
  cell_to_values = lazy_map(Reindex(bgcell_to_values),cell_to_bgcell)
  return cell_to_values
end

function extract_dualized_cell_values(
  ::SubFacetBoundaryTriangulation,
  φh::CellField,
  caches,
)
  @assert isa(DomainStyle(φh),ReferenceDomain)
  face_to_fields = lazy_map(Reindex(CellData.get_data(φh)),caches.face_to_φscell)
  face_to_values = lazy_map(evaluate,face_to_fields,caches.face_to_φrcoords)
  return lazy_map(real,face_to_values)
end

function extract_dualized_cell_values(
  bgmodel::DiscreteModel,
  φh::CellField,
)
  @assert isa(DomainStyle(φh),ReferenceDomain)
  strian = get_triangulation(φh)
  @check get_background_model(strian) === bgmodel
  D = num_cell_dims(bgmodel)
  scell_to_polys = expand_cell_data(get_polytopes(strian),get_cell_type(strian))
  scell_to_rcoords = lazy_map(get_vertex_coordinates,scell_to_polys)
  scell_to_values = lazy_map(evaluate,CellData.get_data(φh),scell_to_rcoords)
  # When φh has complex-valued DOFs, we define the geometry using the real part.
  scell_to_rvalues = lazy_map(real,scell_to_values)
  return extend(scell_to_rvalues,get_glue(strian,Val(D)).mface_to_tface)
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

function FESpaces.get_cell_fe_data(fun,f,ttrian::DifferentiableTriangulationView)
  lazy_map(Reindex(FESpaces.get_cell_fe_data(fun,f,ttrian.parent)),ttrian.cell_to_parent_cell)
end

# AppendedTriangulation
#
# When cutting an embedded domain, we will usually end up with an AppendedTriangulation
# containing
#   a) a regular triangulation with the IN/OUT cells
#   b) a SubFaceTriangulation / SubFacetBoundaryTriangulation with the CUT cells
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

function Geometry.is_change_possible(
  strian::AppendedTriangulation{Dc,Dp},
  ttrian::DifferentiableAppendedTriangulation{Dc,Dp}
) where {Dc,Dp}
  return is_change_possible(strian.a,ttrian.a) && (strian.b === ttrian.b || is_change_possible(strian.b,ttrian.b))
end

function Geometry.best_target(
  strian::AppendedTriangulation{Dc,Dp},
  ttrian::DifferentiableAppendedTriangulation{Dc,Dp}
) where {Dc,Dp}
  return ttrian
end

for tdomain in (:ReferenceDomain,:PhysicalDomain)
  for sdomain in (:ReferenceDomain,:PhysicalDomain)
    @eval begin
      function CellData.change_domain(
        a::CellField,strian::AppendedTriangulation{Dc,Dp},::$sdomain,ttrian::DifferentiableAppendedTriangulation{Dc,Dp},::$tdomain
      ) where {Dc,Dp}
        @assert is_change_possible(strian,ttrian)
        b = change_domain(a,$(tdomain)())
        return CellData.similar_cell_field(a,CellData.get_data(b),ttrian,$(tdomain)())
      end
    end
  end
end

function FESpaces.get_cell_fe_data(fun,f,ttrian::DifferentiableAppendedTriangulation)
  lazy_append(FESpaces.get_cell_fe_data(fun,f,ttrian.a),FESpaces.get_cell_fe_data(fun,f,ttrian.b))
end

const DifferentiableTriangulationTypes = Union{
  DifferentiableTriangulation,
  DifferentiableTriangulationView,
  DifferentiableAppendedTriangulation
}

function FESpaces._change_argument(
  op,f,trian::DifferentiableTriangulationTypes,uh
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
