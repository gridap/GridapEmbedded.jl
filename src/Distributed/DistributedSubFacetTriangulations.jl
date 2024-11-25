
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
