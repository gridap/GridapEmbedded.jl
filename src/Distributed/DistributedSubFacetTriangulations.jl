
const DistributedSubFacetTriangulation{Df,Dc} = DistributedTriangulation{Df,Dc,<:AbstractArray{<:Union{SubFacetTriangulation{Df,Dc},TriangulationView{Df,Dc,<:SubFacetTriangulation{Df,Dc}}}}}

# Each cut facet belongs to the background cell containing it. So we can generate 
# ownership information for the cut facets from the background cell gids.
function GridapDistributed.generate_cell_gids(
  trian::DistributedSubFacetTriangulation{Df,Dc},
) where {Df,Dc}
  model = get_background_model(trian)
  cgids = get_cell_gids(model)

  n_lfacets, bgcell_to_lfacets = map(local_views(trian)) do trian
    model = get_background_model(trian)
    lfacet_to_bgcell = get_glue(trian,Val(Dc)).tface_to_mface
    n_lfacets = length(lfacet_to_bgcell)

    ptrs  = zeros(Int32,num_cells(model)+1)
    for bgcell in lfacet_to_bgcell
      ptrs[bgcell+1] += 1
    end
    Arrays.length_to_ptrs!(ptrs)
    @assert ptrs[end] == n_lfacets+1

    data = zeros(Int32,n_lfacets)
    for (lfacet,bgcell) in enumerate(lfacet_to_bgcell)
      data[ptrs[bgcell]] = lfacet
      ptrs[bgcell] += 1
    end
    Arrays.rewind_ptrs!(ptrs)

    return n_lfacets, Table(data,ptrs)
  end |> tuple_of_arrays

  return GridapDistributed.generate_gids(
    cgids, bgcell_to_lfacets, n_lfacets
  )
end

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
