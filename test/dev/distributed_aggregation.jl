
using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays

using GridapEmbedded: aggregate

function generate_dumbell(ranks, np, nc, ng)
  L = 1
  p0 = Point(0.0,0.0)
  pmin = Point(-L/4,-L/8)
  pmax = Point(L/4,3*L/8)

  R = 0.2
  h = L/(2*nc)

  geo1 = disk(R,x0=p0)
  geo2 = quadrilateral(;x0=Point(-L/4+h/4,5*h/4),
                        d1=VectorValue(8*h+h/2,0.0),
                        d2=VectorValue(0.0,5*h+h/2))
  geo3 = quadrilateral(;x0=Point(-L/4+3*h/4,7*h/4),
                        d1=VectorValue(8*h-h/2,0.0),
                        d2=VectorValue(0.0,5*h-h/2))
  geo = union(geo1,intersect(geo2,!geo3))
  bgmodel = CartesianDiscreteModel(ranks,np,pmin,pmax,(nc,nc);ghost=(ng,ng))
  return bgmodel, geo
end

function exchange_impl!(vector_partition,cache)
  buffer_snd = map(vector_partition,cache) do values, cache
    local_indices_snd = cache.local_indices_snd
    for (p,lid) in enumerate(local_indices_snd.data)
      cache.buffer_snd.data[p] = values[lid]
    end
    cache.buffer_snd
  end
  neighbors_snd, neighbors_rcv, buffer_rcv = map(cache) do cache
    cache.neighbors_snd, cache.neighbors_rcv, cache.buffer_rcv
  end |> tuple_of_arrays
  graph = ExchangeGraph(neighbors_snd,neighbors_rcv)
  t = exchange!(buffer_rcv,buffer_snd,graph)
  return t
end

function find_optimal_roots!(lcell_to_root,lcell_to_path_length,lcell_to_bbox_diam,cell_indices)
  # Bring all root candidates to the owner of the cut cell
  roots_cache = PartitionedArrays.p_vector_cache(lcell_to_root,cell_indices)
  t1 = exchange_impl!(lcell_to_root,roots_cache)
  lengths_cache = PartitionedArrays.p_vector_cache(lcell_to_path_length,cell_indices)
  t2 = exchange_impl!(lcell_to_path_length,lengths_cache)
  bb_diams_cache = PartitionedArrays.p_vector_cache(lcell_to_bbox_diam,cell_indices)
  t3 = exchange_impl!(lcell_to_bbox_diam,bb_diams_cache)
  lcell_to_owner = map(copy∘local_to_owner,cell_indices)
  owners_cache = PartitionedArrays.p_vector_cache(lcell_to_owner,cell_indices)
  wait(t1)
  wait(t2)
  wait(t3)

  # Select the optimal root for each local cut cell
  map(
    lcell_to_root,lcell_to_path_length,lcell_to_bbox_diam,
      lcell_to_owner,roots_cache,lengths_cache,bb_diams_cache
  ) do lcell_to_root,lcell_to_path_length,lcell_to_bbox_diam,
      lcell_to_owner,roots_cache,lengths_cache,bb_diams_cache

    lids_rcv = roots_cache.local_indices_rcv
    lengths_rcv = lengths_cache.buffer_rcv
    bb_diams_rcv = bb_diams_cache.buffer_rcv
    roots_rcv = roots_cache.buffer_rcv

    for (k,nbor) in enumerate(roots_cache.neighbors_rcv)
      for (lcell,root,len,bb_diam) in zip(lids_rcv[k],roots_rcv[k],lengths_rcv[k],bb_diams_rcv[k])
        root == 0 && continue
        if lcell_to_path_length[lcell] > len
          lcell_to_root[lcell] = root
          lcell_to_owner[lcell] = nbor
        elseif lcell_to_path_length[lcell] == len
          if ( lcell_to_bbox_diam[lcell] > bb_diam ) &
              !isapprox(lcell_to_bbox_diam[lcell],bb_diam,atol=1.0e-9)
            lcell_to_root[lcell] = root
            lcell_to_owner[lcell] = nbor
          elseif isapprox(lcell_to_bbox_diam[lcell],bb_diam,atol=1.0e-9) &
              ( lcell_to_root[lcell] > root )
            lcell_to_root[lcell] = root
            lcell_to_owner[lcell] = nbor
          end
        end
      end
    end

  end

  # Scatter the optimal roots and values
  # Technically, we do not need to scatter the values, it can be removed after 
  # we are done debugging
  t1 = consistent!(PVector(lcell_to_owner,cell_indices,owners_cache))
  t2 = consistent!(PVector(lcell_to_root,cell_indices,roots_cache))
  
  wait(t1)

  agg_cell_indices = map(cell_indices,lcell_to_owner) do cell_indices, lcell_to_owner
    LocalIndices(global_length(cell_indices),part_id(cell_indices),local_to_global(cell_indices),lcell_to_owner)
  end

  wait(t2)

  return lcell_to_root, lcell_to_owner, agg_cell_indices
end

distribute = PartitionedArrays.DebugArray
np = (3,1)
ranks = distribute(LinearIndices((prod(np),)))
bgmodel, geo = generate_dumbell(ranks, np, 9, 2)
cutgeo = cut(bgmodel, geo)

cell_indices = partition(get_cell_gids(bgmodel))

strategy = AggregateCutCellsByThreshold(1.0)
lcell_to_root, lcell_to_path_length, lcell_to_bbox_diam =
  map(local_views(cutgeo),cell_indices) do cutgeo,cell_indices
    lid_to_gid = local_to_global(cell_indices)
    aggregate(strategy,cutgeo,geo,lid_to_gid,IN)
  end |> tuple_of_arrays

lcell_to_inconsistent_root = map(copy,lcell_to_root)

lcell_to_root, lcell_to_owner, agg_cell_indices =
  find_optimal_roots!(lcell_to_root,lcell_to_path_length,lcell_to_bbox_diam,cell_indices);

# Creating aggregate-conforming cell partition

gids = get_cell_gids(bgmodel)
ocell_to_root = map(lcell_to_root,own_to_local(gids)) do agg,o_to_l
  map(Reindex(agg),o_to_l)
end

writevtk(EmbeddedBoundary(cutgeo),"data/bnd");
writevtk(
  Triangulation(bgmodel), "data/dumbell_aggregates", 
  celldata = ["aggregate" => ocell_to_root],
);

map(ranks,
    local_views(bgmodel),
    lcell_to_root,
    lcell_to_inconsistent_root,
    lcell_to_owner,
    lcell_to_path_length,
    lcell_to_bbox_diam) do r,bgmodel,
                          lcell_to_root,
                          lcell_to_iroot,
                          lcell_to_owner,
                          lcell_to_path_length,
                          lcell_to_bbox_diam
  writevtk(
    Triangulation(bgmodel), "data/dumbell_aggregates_$(r)", 
    celldata = [ "roots"              => lcell_to_root,
                 "inconsistent roots" => lcell_to_iroot,
                 "owners"             => lcell_to_owner,
                 "length"             => lcell_to_path_length,
                 "diameter"           => lcell_to_bbox_diam ],
  );
end