
const DistributedSubFacetTriangulation{Df,Dc} = DistributedTriangulation{Df,Dc,<:AbstractArray{<:Union{SubFacetTriangulation{Df,Dc},TriangulationView{Df,Dc,<:SubFacetTriangulation{Df,Dc}}}}}

# Each cut facet belongs to the background cell containing it. So we can generate 
# ownership information for the cut facets from the background cell gids.
function GridapDistributed.generate_cell_gids(
  trian::DistributedSubFacetTriangulation{Df,Dc},
) where {Df,Dc}
  model = get_background_model(trian)
  cgids = get_cell_gids(model)
  
  n_lfacets = map(num_cells,local_views(trian))
  first_gid = scan(+,n_lfacets,type=:exclusive,init=one(eltype(n_lfacets)))
  n_facets = reduce(+,n_lfacets,init=zero(eltype(n_lfacets)))

  fgids = map(local_views(trian),partition(cgids),first_gid,n_lfacets) do trian, cgids, first_gid, n_lfacets
    glue = get_glue(trian,Val(Dc)) # Glue from cut facets to background cells 
    facet_to_bgcell = glue.tface_to_mface

    facet_to_gid = collect(first_gid:(first_gid+n_lfacets-1))
    facet_to_owner = local_to_owner(cgids)[facet_to_bgcell]
    LocalIndices(n_facets,part_id(cgids),facet_to_gid,facet_to_owner)
  end
  return PRange(fgids)
end

# function GridapDistributed.generate_cell_gids(
#   trian::DistributedSubFacetTriangulation{Df,Dc},
# ) where {Df,Dc}
#   model = get_background_model(trian)
#   ctrians = map(local_views(trian),local_views(model)) do trian, model
#     glue = get_glue(trian,Val(Dc)) # Glue from cut facets to background cells 
#     facet_to_bgcell = glue.tface_to_mface
#     Triangulation(model,facet_to_bgcell)
#   end
#   ctrian = GridapDistributed.DistributedTriangulation(ctrians,model)
#   return GridapDistributed.generate_cell_gids(ctrian)
# end

function GridapDistributed.add_ghost_cells(
  trian::DistributedSubFacetTriangulation{Df,Dc},
) where {Df,Dc}

  # In this case, we already have all ghost facets
  if eltype(local_views(trian)) <: SubFacetTriangulation
    return trian
  end

  # First, we create a new Triangulation containing all the cut facets
  model = get_background_model(trian)
  bgtrians, facet_to_bgfacet = map(local_views(trian)) do trian
    @assert isa(trian,TriangulationView)
    trian.parent, trian.cell_to_parent_cell
  end |> tuple_of_arrays
  bgtrian = DistributedTriangulation(bgtrians,model)
  fgids = partition(generate_cell_gids(bgtrian))

  # Exchange info about cut facets
  inside_facets = map(fgids,facet_to_bgfacet) do fgids, facet_to_bgfacet
    inside_facets = falses(local_length(fgids))
    inside_facets[facet_to_bgfacet] .= true
    return inside_facets
  end
  wait(consistent!(PVector(inside_facets,fgids))) # Exchange information

  # Return ghosted Triangulation
  covers_all = reduce(&,map(all,inside_facets),init=true)
  if covers_all
    ghosted_trian = bgtrian
  else
    ghosted_trian = DistributedTriangulation(
      map(TriangulationView,bgtrians,inside_facets), model
    )
  end
  return ghosted_trian
end

function num_cells(trian::DistributedSubFacetTriangulation)
  model = get_background_model(trian)
  Dc = num_cell_dims(model)
  gids = get_face_gids(model,Dc)
  n_loc_ocells = map(local_views(trian),partition(gids)) do trian, gids
    glue = get_glue(trian,Val(Dc))
    @assert isa(glue,FaceToFaceGlue)
    tcell_to_mcell = glue.tface_to_mface
    if isa(tcell_to_mcell,IdentityVector)
      own_length(gids)
    else
      mcell_to_owned = local_to_own(gids)
      is_owned(mcell) = !iszero(mcell_to_owned[mcell])
      sum(is_owned,tcell_to_mcell;init=0)
    end
  end
  return sum(n_loc_ocells)
end
