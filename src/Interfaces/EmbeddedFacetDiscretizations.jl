
"""
    struct EmbeddedFacetDiscretization{Dc,Dp,T} <: AbstractEmbeddedDiscretization

This structure contains all the required information to build integration `Triangulations` 
for a cut model boundary.

## Constructors

    cut_facets(cutter::Cutter,background,geom)

## Properties

- `bgmodel::DiscreteModel`: the background mesh
- `geo::CSG.Geometry`: the geometry used to cut the background mesh
- `subfacets::SubFacetData`: collection of cut facets, attached to the background mesh
- `ls_to_facet_to_inoutcut::Vector{Vector{Int8}}`: list of IN/OUT/CUT states for each facet 
   in the background mesh, for each node in the geometry tree.
- `ls_to_subfacet_to_inoutcut::Vector{Vector{Int8}}`: list of IN/OUT/CUT states for each subfacet 
   in the cut part of the mesh, for each node in the geometry tree.

## Methods

- [`BoundaryTriangulation(cut::EmbeddedFacetDiscretization,in_or_out)`](@ref)
- [`SkeletonTriangulation(cut::EmbeddedFacetDiscretization,in_or_out)`](@ref)

"""
struct EmbeddedFacetDiscretization{Dc,Dp,T} <: AbstractEmbeddedDiscretization
  bgmodel::DiscreteModel{Dp,Dp}
  ls_to_facet_to_inoutcut::Vector{Vector{Int8}}
  subfacets::SubCellData{Dc,Dp,T}
  ls_to_subfacet_to_inout::Vector{Vector{Int8}}
  oid_to_ls::Dict{UInt,Int}
  geo::CSG.Geometry
end

function get_background_model(cut::EmbeddedFacetDiscretization)
  cut.bgmodel
end

function get_geometry(cut::EmbeddedFacetDiscretization)
  cut.geo
end

function SkeletonTriangulation(cut::EmbeddedFacetDiscretization)
  SkeletonTriangulation(cut,PHYSICAL_IN)
end

"""
    SkeletonTriangulation(cut::EmbeddedFacetDiscretization[, in_or_out=PHYSICAL_IN])
"""
function SkeletonTriangulation(
  cut::EmbeddedFacetDiscretization,
  in_or_out)

  facets = SkeletonTriangulation(cut.bgmodel)
  geo = cut.geo
  SkeletonTriangulation(facets,cut,in_or_out,geo)
end

function SkeletonTriangulation(cut::EmbeddedFacetDiscretization,name::String)
  SkeletonTriangulation(cut,PHYSICAL_IN,name)
end

function SkeletonTriangulation(cut::EmbeddedFacetDiscretization,geo::CSG.Geometry)
  SkeletonTriangulation(cut,PHYSICAL_IN,geo)
end

function SkeletonTriangulation(
  cut::EmbeddedFacetDiscretization,
  in_or_out,
  name::String)

  facets = SkeletonTriangulation(cut.bgmodel)
  geo = get_geometry(cut.geo,name)
  SkeletonTriangulation(facets,cut,in_or_out,geo)
end

function SkeletonTriangulation(
  cut::EmbeddedFacetDiscretization,
  in_or_out,
  geo::CSG.Geometry)

  facets = SkeletonTriangulation(cut.bgmodel)
  SkeletonTriangulation(facets,cut,in_or_out,geo)
end

function SkeletonTriangulation(
  facets::SkeletonTriangulation,
  cut::EmbeddedFacetDiscretization,
  in_or_out,
  geo::CSG.Geometry)

  facets1 = facets.⁺
  facets2 = facets.⁻
  trian1 = BoundaryTriangulation(facets1,cut,in_or_out,geo)
  trian2 = BoundaryTriangulation(facets2,cut,in_or_out,geo)
  SkeletonTriangulation(trian1,trian2)
end

function BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization;
  tags=nothing)
  BoundaryTriangulation(cut,PHYSICAL_IN;tags=tags)
end

"""
    BoundaryTriangulation(cut::EmbeddedFacetDiscretization[, in_or_out=PHYSICAL_IN; tags=nothing])
"""
function BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization,
  in_or_out;
  tags=nothing)

  facets = BoundaryTriangulation(cut.bgmodel;tags=tags)
  geo = cut.geo
  BoundaryTriangulation(facets,cut,in_or_out,geo)
end

function BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization,
  name::String;
  tags=nothing)
  BoundaryTriangulation(cut,PHYSICAL_IN,name;tags=tags)
end

function BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization,
  geo::CSG.Geometry;
  tags=nothing)
  BoundaryTriangulation(cut,PHYSICAL_IN,geo;tags=tags)
end

function BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization,
  in_or_out,
  name::String;
  tags=nothing)

  facets = BoundaryTriangulation(cut.bgmodel;tags=tags)
  geo = get_geometry(cut.geo,name)
  BoundaryTriangulation(facets,cut,in_or_out,geo)
end

function BoundaryTriangulation(
  cut::EmbeddedFacetDiscretization,
  in_or_out,
  geo::CSG.Geometry;
  tags=nothing)

  facets = BoundaryTriangulation(cut.bgmodel;tags=tags)
  BoundaryTriangulation(facets,cut,in_or_out,geo)
end

function BoundaryTriangulation(
  facets::BoundaryTriangulation,
  cut::EmbeddedFacetDiscretization,
  in_or_out::Tuple,
  geo::CSG.Geometry)

  a = BoundaryTriangulation(facets,cut,in_or_out[1],geo)
  b = BoundaryTriangulation(facets,cut,in_or_out[2],geo)
  iszero(num_cells(a)) ? b : lazy_append(a,b)
end

function BoundaryTriangulation(
  facets::BoundaryTriangulation,
  cut::EmbeddedFacetDiscretization,
  in_or_out::Integer,
  geo::CSG.Geometry)

  bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(cut,geo)
  bgfacet_to_mask = lazy_map(isequal(in_or_out), bgfacet_to_inoutcut)
  _restrict_boundary_triangulation(cut.bgmodel,facets,bgfacet_to_mask)
end

function BoundaryTriangulation(
  facets::BoundaryTriangulation,
  cut::EmbeddedFacetDiscretization,
  in_or_out::ActiveInOrOut,
  geo::CSG.Geometry)

  bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(cut,geo)
  bgfacet_to_mask = lazy_map( a->a==CUT || a==in_or_out.in_or_out, bgfacet_to_inoutcut)
  _restrict_boundary_triangulation(cut.bgmodel,facets,bgfacet_to_mask)
end

function BoundaryTriangulation(
  _facets::BoundaryTriangulation,
  cut::EmbeddedFacetDiscretization,
  in_or_out::CutInOrOut,
  geo::CSG.Geometry)

  bgfacet_to_inoutcut = compute_bgfacet_to_inoutcut(cut,geo)
  bgfacet_to_mask = lazy_map(isequal(CUT), bgfacet_to_inoutcut)
  facets = _restrict_boundary_triangulation(cut.bgmodel,_facets,bgfacet_to_mask)

  facet_to_bgfacet = facets.glue.face_to_bgface
  n_bgfacets = num_facets(cut.bgmodel)
  bgfacet_to_facet = zeros(Int,n_bgfacets)
  bgfacet_to_facet[facet_to_bgfacet] .= 1:length(facet_to_bgfacet)

  subfacet_to_inoutcut = lazy_map(Reindex(bgfacet_to_inoutcut),cut.subfacets.cell_to_bgcell)
  _subfacet_to_facet = lazy_map(Reindex(bgfacet_to_facet),cut.subfacets.cell_to_bgcell)

  subfacet_to_inout = compute_subfacet_to_inout(cut,geo)
  pred(a,b,c) = !iszero(c) && a==CUT && b==in_or_out.in_or_out
  mask = lazy_map( pred, subfacet_to_inoutcut, subfacet_to_inout, _subfacet_to_facet )
  newsubfacets = findall(mask)
  subfacets = SubCellData(cut.subfacets,newsubfacets)
  subfacet_to_facet = bgfacet_to_facet[subfacets.cell_to_bgcell]

  SubFacetBoundaryTriangulation(facets,subfacets,subfacet_to_facet)
end

function _restrict_boundary_triangulation(model,facets,bgfacet_to_mask)
  facet_to_bgfacet = facets.glue.face_to_bgface
  facet_to_mask = lazy_map(Reindex(bgfacet_to_mask),facet_to_bgfacet)
  n_bgfacets = length(bgfacet_to_mask)
  bgfacet_to_mask2 = fill(false,n_bgfacets)
  bgfacet_to_mask2[facet_to_bgfacet] .= facet_to_mask

  BoundaryTriangulation(model,bgfacet_to_mask2,facets.glue.bgface_to_lcell)
end

function compute_bgfacet_to_inoutcut(cut::EmbeddedFacetDiscretization,geo::CSG.Geometry)

  tree = get_tree(geo)

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    ls = cut.oid_to_ls[oid]
    cell_to_inoutcut = cut.ls_to_facet_to_inoutcut[ls]
    cell_to_inoutcut, name, meta
  end

  newtree = replace_data(identity,conversion,tree)
  compute_inoutcut(newtree)
end

function compute_subfacet_to_inout(cut::EmbeddedFacetDiscretization,geo::CSG.Geometry)

  tree = get_tree(geo)

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    ls = cut.oid_to_ls[oid]
    cell_to_inoutcut = cut.ls_to_subfacet_to_inout[ls]
    cell_to_inoutcut, name, meta
  end

  newtree = replace_data(identity,conversion,tree)
  compute_inoutcut(newtree)
end

"""
    struct SubFacetBoundaryTriangulation{Dc,Dp,T} <: Triangulation{Dc,Dp}

Triangulation of cut facets from the background mesh, i.e each of the facets 
in this triangulation is part of a background facet that has been cut by the geometry.

This differs from the the `SubFacetTriangulation` in that the facets in the `SubFacetTriangulation` 
are not cut background facets, but rather subfacets on the interior of a background cell.

They result from calling `Boundary` or `Skeleton` on an `EmbeddedFacetDiscretization` object.

    BoundaryTriangulation(cut::EmbeddedFacetDiscretization,in_or_out;tags=nothing)
    SkeletonTriangulation(cut::EmbeddedFacetDiscretization,in_or_out)

"""
struct SubFacetBoundaryTriangulation{Dc,Dp,T} <: Triangulation{Dc,Dp}
  facets::BoundaryTriangulation{Dc,Dp}
  subfacets::SubCellData{Dc,Dp,T}
  subfacet_to_facet::AbstractArray
  reffes::Vector{LagrangianRefFE{Dc}}
  cell_types::Vector{Int8}
  cell_ids
  cell_normals
  cell_ref_map
  subgrid

  function SubFacetBoundaryTriangulation(
    facets::BoundaryTriangulation{Dc,Dp},
    subfacets::SubCellData{Dc,Dp,T},
    subfacet_to_facet::AbstractArray) where {Dc,Dp,T}

    reffe = LagrangianRefFE(Float64,Simplex(Val{Dc}()),1)
    cell_types = fill(Int8(1),length(subfacets.cell_to_points))
    reffes = [reffe]
    glue = get_glue(facets,Val(Dc+1))
    subgrid = UnstructuredGrid(subfacets)

    cell_ids = lazy_map(Reindex(glue.tface_to_mface),subfacet_to_facet)
    cell_normals = lazy_map(Reindex(get_facet_normal(facets)),subfacet_to_facet)
    subfacet_to_facet_map = _setup_cell_ref_map(subfacets,subgrid)
    face_ref_map = lazy_map(Reindex(glue.tface_to_mface_map),subfacet_to_facet)
    cell_ref_map = lazy_map(∘,face_ref_map,subfacet_to_facet_map)

    new{Dc,Dp,T}(
      facets,
      subfacets,
      subfacet_to_facet,
      reffes,
      cell_types,
      cell_ids,
      cell_normals,
      cell_ref_map,
      subgrid)
  end
end

function get_background_model(a::SubFacetBoundaryTriangulation)
  get_background_model(a.facets)
end

function get_active_model(a::SubFacetBoundaryTriangulation)
  msg = """
  This is not implemented, but also not needed in practice.
  Embedded Grids implemented for integration, not interpolation.
  """
  @notimplemented  msg
end

function get_grid(a::SubFacetBoundaryTriangulation)
  a.subgrid
end

function get_glue(a::SubFacetBoundaryTriangulation{Dc},::Val{D}) where {Dc,D}
  if D == Dc
    tface_to_mface = a.subfacets.cell_to_bgcell
    tface_to_mface_map = _setup_cell_ref_map(a.subfacets,a.subgrid)
    FaceToFaceGlue(tface_to_mface,tface_to_mface_map,nothing)
  elseif D-1 == Dc
    tface_to_mface = a.cell_ids
    tface_to_mface_map = a.cell_ref_map
    FaceToFaceGlue(tface_to_mface,tface_to_mface_map,nothing)
  else
    nothing
  end
end

function get_facet_normal(trian::SubFacetBoundaryTriangulation)
  trian.cell_normals
end

