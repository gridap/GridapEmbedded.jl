
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
  x = VectorValue(1.0,0.0)
  h = L/(2*nc)

  geo1 = disk(R,x0=p0)
  geo2 = quadrilateral(;x0=Point(-L/4+h/4,5*h/4),
                        d1=VectorValue(6*h+h/2,0.0),
                        d2=VectorValue(0.0,4*h+h/2))
  geo3 = quadrilateral(;x0=Point(-L/4+3*h/4,7*h/4),
                        d1=VectorValue(6*h-h/2,0.0),
                        d2=VectorValue(0.0,4*h-h/2))
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

function find_optimal_root!(lcell_to_root,lcell_to_value,cell_indices)
  # Bring all root candidates to the owner of the cut cell
  roots_cache = PartitionedArrays.p_vector_cache(lcell_to_lroot,cell_indices)
  t1 = exchange_impl!(lcell_to_lroot,roots_cache)
  values_cache = PartitionedArrays.p_vector_cache(lcell_to_value,cell_indices)
  t2 = exchange_impl!(lcell_to_value,values_cache)
  wait(t1)
  wait(t2)

  # Select the optimal root for each local cut cell
  map(lcell_to_root,lcell_to_value,roots_cache,values_cache) do lcell_to_root, lcell_to_value, roots_cache, values_cache
    local_indices_rcv = roots_cache.local_indices_rcv.data
    values_rcv = values_cache.buffer_rcv.data
    roots_rcv = roots_cache.buffer_rcv.data
    for (p, lcell) in enumerate(local_indices_rcv)
      value, root = values_rcv[p], roots_rcv[p]
      if lcell_to_value[lcell] > value # Take the minimum
        lcell_to_root[lcell] = root
        lcell_to_value[lcell] = value
      end
    end
  end

  # Scatter the optimal roots and values
  # Technically, we do not need to scatter the values, it can be removed after 
  # we are done debugging
  t1 = consistent!(PVector(lcell_to_root,cell_indices,roots_cache))
  t2 = consistent!(PVector(lcell_to_value,cell_indices,values_cache))
  wait(t1)
  wait(t2)

  return lcell_to_root, lcell_to_value
end


distribute = PartitionedArrays.DebugArray
np = (2,1)
ranks = distribute(LinearIndices((prod(np),)))
bgmodel, geo = generate_dumbell(ranks, np, 8, 3)
cutgeo = cut(bgmodel, geo)

cell_indices = partition(get_cell_gids(bgmodel))

strategy = AggregateCutCellsByThreshold(1.0)
lcell_to_lroot = map(local_views(cutgeo)) do cutgeo
  aggregate(strategy,cutgeo,geo,IN)
end

# This would ideally come from the serial algorithm, but we simulate it here
lcell_to_value = map(local_views(bgmodel),lcell_to_lroot) do bgmodel, lcell_to_lroot
  lcell_centroid = lazy_map(mean, get_cell_coordinates(bgmodel))
  lcell_to_value = fill(Inf, num_cells(bgmodel))
  for (lcell, lroot) in enumerate(lcell_to_lroot)
    if !iszero(lroot)
      value = norm(lcell_centroid[lcell] - lcell_centroid[lroot])
      lcell_to_value[lcell] = value
    end
  end
  return lcell_to_value
end
lcell_to_value_copy = map(copy,lcell_to_value)

# This can also be done in place, but we duplicate infor for now
lcell_to_root = map(lcell_to_lroot,cell_indices) do lcell_to_lroot, cell_indices
  lcell_to_root = fill(0, length(cell_indices))
  lid_to_gid = local_to_global(cell_indices)
  for i in eachindex(lcell_to_root)
    if !iszero(lcell_to_lroot[i])
      lcell_to_root[i] = lid_to_gid[lcell_to_lroot[i]]
    end
  end
  return lcell_to_root
end

lcell_to_root, lcell_to_value = find_optimal_root!(lcell_to_root,lcell_to_value,cell_indices);

display(lcell_to_root)

gids = get_cell_gids(bgmodel)
ocell_to_root = map(lcell_to_root,own_to_local(gids)) do agg,o_to_l
  map(Reindex(agg),o_to_l)
end
ocell_to_lroot = map(lcell_to_lroot,own_to_local(gids)) do agg,o_to_l
  map(Reindex(agg),o_to_l)
end

writevtk(EmbeddedBoundary(cutgeo),"data/bnd");
writevtk(
  Triangulation(bgmodel), "data/dumbell_aggregates", 
  celldata = ["aggregate" => ocell_to_root, "local_aggregates" => ocell_to_lroot],
);

map(ranks,local_views(bgmodel),lcell_to_lroot,lcell_to_value_copy) do r,bgmodel, lcell_to_lroot, lcell_to_value
  writevtk(
    Triangulation(bgmodel), "data/dumbell_aggregates_$(r)", 
    celldata = ["values" => lcell_to_value, "roots" => lcell_to_lroot],
  );
end
