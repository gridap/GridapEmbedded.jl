
struct DistributedEmbeddedDiscretization{A,B} <: GridapType
  discretizations::A
  model::B
  function DistributedEmbeddedDiscretization(
    discretizations::AbstractArray{<:AbstractEmbeddedDiscretization},
    model::DistributedDiscreteModel
  )
    A = typeof(discretizations)
    B = typeof(model)
    new{A,B}(discretizations,model)
  end
end

local_views(a::DistributedEmbeddedDiscretization) = a.discretizations

get_background_model(a::DistributedEmbeddedDiscretization) = a.model

function get_geometry(a::DistributedEmbeddedDiscretization)
  geometries = map(get_geometry,local_views(a))
  distributed_geometry(geometries)
end

# Needed for dispatching between analytical geometries and discrete geometries
function distributed_geometry(geometries::AbstractArray{<:CSG.Geometry})
  PartitionedArrays.getany(geometries)
end

function cut(bgmodel::DistributedDiscreteModel,args...)
  cut(LevelSetCutter(),bgmodel,args...)
end

function cut(cutter::Cutter,bgmodel::DistributedDiscreteModel{Dc},args...) where Dc
  gids = get_face_gids(bgmodel,Dc)
  cuts = map(local_views(bgmodel)) do bgmodel
    cut(cutter,bgmodel,args...)
  end
  @notimplementedif !isconsistent_bgcell_to_inoutcut(cuts,partition(gids))
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

function cut_facets(bgmodel::DistributedDiscreteModel,args...)
  cut_facets(LevelSetCutter(),bgmodel,args...)
end

function cut_facets(cutter::Cutter,bgmodel::DistributedDiscreteModel{Dc},args...) where Dc
  gids = get_face_gids(bgmodel,Dc-1)
  cuts = map(local_views(bgmodel)) do bgmodel
    cut_facets(cutter,bgmodel,args...)
  end
  @notimplementedif !isconsistent_bgcell_to_inoutcut(cuts,partition(gids))
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

# Note on distributed triangulations:
#
# - We allow for more one argument, `portion`, which allows the user to filter
# some of the cells/faces. In particular, this is used to remove ghosts from the
# local triangulations.
# - The default for `portion` is `NoGhost()`, wich filters out all ghost cells, except
# when we have the argument `in_or_out`.

function Triangulation(
  cutgeo::DistributedEmbeddedDiscretization,in_or_out::ActiveInOrOut,args...
)
  Triangulation(WithGhost(),cutgeo,in_or_out,args...)
end

for TT in (:Triangulation,:SkeletonTriangulation,:BoundaryTriangulation,:EmbeddedBoundary,:GhostSkeleton)
  @eval begin
    function $TT(cutgeo::DistributedEmbeddedDiscretization,args...)
      $TT(NoGhost(),cutgeo,args...)
    end

    function $TT(portion,cutgeo::DistributedEmbeddedDiscretization,args...)
      model = get_background_model(cutgeo)
      gids  = get_cell_gids(model)
      trians = map(local_views(cutgeo),partition(gids)) do cutgeo, gids
        $TT(portion,gids,cutgeo,args...)
      end
      DistributedTriangulation(trians,model)
    end

    function $TT(portion,gids::AbstractLocalIndices,cutgeo::AbstractEmbeddedDiscretization,args...)
      trian = $TT(cutgeo,args...)
      filter_cells_when_needed(portion,gids,trian)
    end
  end
end

# TODO: This should go to GridapDistributed
function remove_ghost_cells(trian::AppendedTriangulation,gids)
  a = remove_ghost_cells(trian.a,gids)
  b = remove_ghost_cells(trian.b,gids)
  iszero(num_cells(a)) && return b
  iszero(num_cells(b)) && return a
  return lazy_append(a,b)
end

function remove_ghost_cells(trian::SubFacetTriangulation{Df,Dc},gids) where {Df,Dc}
  glue  = get_glue(trian,Val{Dc}())
  remove_ghost_cells(glue,trian,gids)
end

function remove_ghost_subfacets(cut::EmbeddedFacetDiscretization,facet_gids)
  bgfacet_mask = map(!iszero,local_to_owner(facet_gids))
  subfacet_mask = map(Reindex(bgfacet_mask),cut.subfacets.cell_to_bgcell)
  new_subfacets = findall(subfacet_mask)
  subfacets = SubCellData(cut.subfacets,new_subfacets)
  ls_to_subfacet_to_inout = map(cut.ls_to_subfacet_to_inout) do sf_to_io
    map(Reindex(sf_to_io),new_subfacets)
  end
  EmbeddedFacetDiscretization(
    cut.bgmodel,
    cut.ls_to_facet_to_inoutcut,
    subfacets,
    ls_to_subfacet_to_inout,
    cut.oid_to_ls,
    cut.geo
  )
end

# Distributed InOutCut flag methods

#     isconsistent_bgcell_to_inoutcut(cut::DistributedEmbeddedDiscretization)
#     isconsistent_bgcell_to_inoutcut(cuts::AbstractArray{<:AbstractEmbeddedDiscretization},indices)
#
# Returns true if the local `ls_to_bgcell_to_inoutcut` arrays are consistent
# accross processors.
function isconsistent_bgcell_to_inoutcut(
  cut::DistributedEmbeddedDiscretization{Dc}
) where Dc
  model = get_background_model(cut)
  gids = get_face_gids(model,Dc)
  isconsistent_bgcell_to_inoutcut(local_views(cut),partition(gids))
end

function isconsistent_bgcell_to_inoutcut(
  cuts::AbstractArray{<:AbstractEmbeddedDiscretization},indices::AbstractArray
)
  get_inoutcut(cut::EmbeddedDiscretization) = Tuple(cut.ls_to_bgcell_to_inoutcut)
  get_inoutcut(cut::EmbeddedFacetDiscretization) = Tuple(cut.ls_to_facet_to_inoutcut)
  ls_to_bgcell_to_inoutcut = tuple_of_arrays(map(get_inoutcut,cuts))
  return isconsistent_bgcell_to_inoutcut(ls_to_bgcell_to_inoutcut,indices)
end

function isconsistent_bgcell_to_inoutcut(
  ls_to_bgcell_to_inoutcut::NTuple{N,<:AbstractArray{<:Vector}},indices::AbstractArray
) where N
  for bgcell_to_inoutcut in ls_to_bgcell_to_inoutcut
    if !isconsistent_bgcell_to_inoutcut(bgcell_to_inoutcut,indices)
      return false
    end
  end
  return true
end

function isconsistent_bgcell_to_inoutcut(
  bgcell_to_inoutcut::AbstractArray{<:Vector},indices::AbstractArray
)
  # TODO: Some allocations can be avoided by going to the low-level communication API
  ref = map(copy,bgcell_to_inoutcut)
  wait(consistent!(PVector(ref,indices)))
  is_consistent = map(bgcell_to_inoutcut,ref) do bgcell_to_inoutcut,ref
    bgcell_to_inoutcut == ref
  end
  return reduce(&,is_consistent,init=true)
end

# TODO: Should we check for consistency here?
function compute_bgfacet_to_inoutcut(bgmodel::DistributedDiscreteModel,args...)
  cutter = LevelSetCutter()
  compute_bgfacet_to_inoutcut(cutter,bgmodel,args...)
end

function compute_bgfacet_to_inoutcut(cutter::Cutter,bgmodel::DistributedDiscreteModel,args...)
  map(local_views(bgmodel)) do bgmodel
    compute_bgfacet_to_inoutcut(cutter,bgmodel,args...)
  end
end

function compute_bgcell_to_inoutcut(cutgeo::DistributedEmbeddedDiscretization,args...)
  map(local_views(cutgeo)) do cutgeo
    compute_bgcell_to_inoutcut(cutgeo,args...)
  end
end

function compute_bgfacet_to_inoutcut(cutgeo::DistributedEmbeddedDiscretization,args...)
  map(local_views(cutgeo)) do cutgeo
    compute_bgfacet_to_inoutcut(cutgeo,args...)
  end
end

# AMR

function compute_redistribute_wights(
  cut::DistributedEmbeddedDiscretization,
  args...)

  geo = get_geometry(cut)
  compute_redistribute_wights(cut,geo,args...)
end

function compute_redistribute_wights(
  cut::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  args...)

  compute_redistribute_wights(compute_bgcell_to_inoutcut(cut,geo),args...)
end

function compute_redistribute_wights(cell_to_inoutcut,in_or_out=IN)
  map(cell_to_inoutcut) do cell_to_inoutcut
    map(cell_to_inoutcut) do inoutcut
      Int( inoutcut ∈ (CUT,in_or_out) )
    end
  end
end

function compute_adaptive_flags(
  cut::DistributedEmbeddedDiscretization,
  args...)

  geo = get_geometry(cut)
  compute_adaptive_flags(cut,geo,args...)
end

function compute_adaptive_flags(
  cut::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  args...)

  compute_adaptive_flags(compute_bgcell_to_inoutcut(cut,geo),args...)
end

function compute_adaptive_flags(cell_to_inoutcut)
  map(cell_to_inoutcut) do c_to_ioc
    flags = zeros(Cint,length(c_to_ioc))
    flags .= nothing_flag
    for (c,ioc) in enumerate(c_to_ioc)
      if ioc == CUT
        flags[c] = refine_flag
      end
    end
    flags
  end
end
