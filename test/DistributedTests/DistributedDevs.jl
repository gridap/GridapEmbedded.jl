module DistributedDevs

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test
using FillArrays


using GridapEmbedded.CSG
using GridapEmbedded.Interfaces: CUT
using Gridap.Geometry
using Gridap.ReferenceFEs

using GridapEmbedded.AgFEM: _touch_aggregated_cells!

using GridapEmbedded.Distributed: DistributedEmbeddedDiscretization

function _aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out=IN)

  bgmodel = get_background_model(cut)
  facet_to_inoutcut = compute_bgfacet_to_inoutcut(bgmodel,geo)
  _distributed_aggregate_by_threshold(strategy.threshold,cut,geo,in_or_out,facet_to_inoutcut)
end


function _distributed_aggregate_by_threshold(threshold,cutgeo,geo,loc,facet_to_inoutcut)
  @assert loc in (IN,OUT)

  cutinorout = loc == IN ? (CUT_IN,IN) : (CUT_OUT,OUT)
  trian = Triangulation(cutgeo,cutinorout,geo)
  model = get_background_model(cutgeo)
  bgtrian = get_triangulation(model)
  cell_to_cut_meas = map(get_cell_measure,local_views(trian),local_views(bgtrian))
  cell_to_meas = map(get_cell_measure,local_views(bgtrian))
  cell_to_unit_cut_meas = map(cell_to_cut_meas,cell_to_meas) do c_to_cm,c_to_m
    lazy_map(/,c_to_cm,c_to_m)
  end

  cell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geo)

  cell_to_coords = map(get_cell_coordinates,local_views(model))
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  cell_to_faces = map(t->get_faces(t,D,D-1),local_views(topo))
  face_to_cells = map(t->get_faces(t,D-1,D),local_views(topo))
  gids = get_cell_gids(model)

  _distributed_aggregate_by_threshold_barrier(
    threshold,cell_to_unit_cut_meas,facet_to_inoutcut,cell_to_inoutcut,
    loc,cell_to_coords,cell_to_faces,face_to_cells,gids)
end


function _distributed_aggregate_by_threshold_barrier(
  threshold,cell_to_unit_cut_meas,facet_to_inoutcut,cell_to_inoutcut,
  loc,cell_to_coords,cell_to_faces,face_to_cells,gids)


  ocell_to_touched = map(cell_to_unit_cut_meas) do c_to_m
    map(≥,c_to_m,Fill(threshold,length(c_to_m)))
  end
  cell_to_touched = _add_ghost_values(ocell_to_touched,gids)

  cell_to_root_centroid = map(cell_to_coords) do cell_to_coords
    map(i->sum(i)/length(i),cell_to_coords)
  end
  PVector(cell_to_root_centroid,partition(gids)) |> consistent! |> wait

  n_cells = map(length,cell_to_touched)
  touched_cells = map(findall,cell_to_touched)

  cell_to_cellin = map(n->zeros(Int32,n),n_cells)
  map(cell_to_cellin,touched_cells,local_to_global(gids)) do c_to_ci,cells,l_to_g
    gcells = lazy_map(Reindex(l_to_g),cells)
    c_to_ci[cells] = gcells
  end

  cell_to_neig = map(n->zeros(Int32,n),n_cells)

  cell_to_root_part = map(collect,local_to_owner(gids))

  c1 = map(array_cache,cell_to_faces)
  c2 = map(array_cache,face_to_cells)

  max_iters = 1
  for iter in 1:max_iters
    all_aggregated = _aggregate_one_step!(c1,c2,gids,
      cell_to_inoutcut,
      cell_to_touched,
      cell_to_neig,
      cell_to_cellin,
      cell_to_root_centroid,
      cell_to_root_part,
      cell_to_faces,
      face_to_cells,
      facet_to_inoutcut,
      loc)

    PVector(cell_to_touched,partition(gids)) |> consistent! |> wait
    PVector(cell_to_neig,partition(gids)) |> consistent! |> wait
    PVector(cell_to_cellin,partition(gids)) |> consistent! |> wait
    PVector(cell_to_root_centroid,partition(gids)) |> consistent! |> wait
    PVector(cell_to_root_part,partition(gids)) |> consistent! |> wait

    reduction!(&,all_aggregated,all_aggregated)
    emit!(all_aggregated,all_aggregated)

    if PartitionedArrays.getany(all_aggregated)
      break
    end
  end

  cell_to_cellin, cell_to_root_part, cell_to_neig
end

function _aggregate_one_step!(c1,c2,gids::PRange,
  cell_to_inoutcut,
  cell_to_touched,
  cell_to_neig,
  cell_to_cellin,
  cell_to_root_centroid,
  cell_to_root_part,
  cell_to_faces,
  face_to_cells,
  facet_to_inoutcut,
  loc)

  map(c1,c2,own_to_local(gids),
    cell_to_inoutcut,
    cell_to_touched,
    cell_to_neig,
    cell_to_cellin,
    cell_to_root_centroid,
    cell_to_root_part,
    local_to_global(gids),
    cell_to_faces,
    face_to_cells,
    facet_to_inoutcut) do c1,c2,own_cells,
        cell_to_inoutcut,
        cell_to_touched,
        cell_to_neig,
        cell_to_cellin,
        cell_to_root_centroid,
        cell_to_root_part,
        cell_to_gcell,
        cell_to_faces,
        face_to_cells,
        facet_to_inoutcut

    _aggregate_one_step!(
      c1,c2,own_cells,
      cell_to_inoutcut,
      cell_to_touched,
      cell_to_neig,
      cell_to_cellin,
      cell_to_root_centroid,
      cell_to_root_part,
      cell_to_gcell,
      cell_to_faces,
      face_to_cells,
      facet_to_inoutcut,
      loc)
  end
end

function _aggregate_one_step!(
  c1,c2,own_cells,
  cell_to_inoutcut,
  cell_to_touched,
  cell_to_neig,
  cell_to_cellin,
  cell_to_root_centroid,
  cell_to_root_part,
  cell_to_gcell,
  cell_to_faces,
  face_to_cells,
  facet_to_inoutcut,
  loc)

  all_aggregated = true
  for cell in own_cells
    if ! cell_to_touched[cell] && cell_to_inoutcut[cell] == CUT
      neigh_cell = _find_best_neighbor_from_centroid_distance(
        c1,c2,cell,
        cell_to_faces,
        face_to_cells,
        cell_to_touched,
        cell_to_root_centroid,
        facet_to_inoutcut,
        loc)
      if neigh_cell > 0
        cellin = cell_to_cellin[neigh_cell]
        centroid = cell_to_root_centroid[neigh_cell]
        part = cell_to_root_part[neigh_cell]
        neigh_gcell = cell_to_gcell[neigh_cell]

        cell_to_neig[cell] = neigh_gcell
        cell_to_cellin[cell] = cellin
        cell_to_root_centroid[cell] = centroid
        cell_to_root_part[cell] = part
      else
        all_aggregated = false
      end
    end
  end
  _touch_aggregated_cells!(cell_to_touched,cell_to_cellin)
  all_aggregated
end

function _find_best_neighbor_from_centroid_distance(
  c1,c2,cell,
  cell_to_faces,
  face_to_cells,
  cell_to_touched,
  cell_to_root_centroid,
  facet_to_inoutcut,
  loc)

  faces = getindex!(c1,cell_to_faces,cell)
  dmin = Inf
  T = eltype(eltype(face_to_cells))
  best_neigh_cell = zero(T)
  for face in faces
    inoutcut = facet_to_inoutcut[face]
    if  inoutcut != CUT && inoutcut != loc
      continue
    end
    neigh_cells = getindex!(c2,face_to_cells,face)
    for neigh_cell in neigh_cells
      if neigh_cell != cell && cell_to_touched[neigh_cell]
        p = cell_to_root_centroid[neigh_cell]
        q = cell_to_root_centroid[cell]
        d = norm(p-q)
        if (1.0+1.0e-9)*d < dmin
          dmin = d
          best_neigh_cell = neigh_cell
        end
      end
    end
  end
  best_neigh_cell
end

function _add_ghost_values(own_v,gids::PRange)
  lens = map(length,local_views(gids))
  eltypes = map(eltype,own_v)
  local_v = map(zeros,eltypes,lens)
  map(local_v,own_v,own_to_local(gids)) do l,o,o_to_l
    l[o_to_l] = o
  end
  PVector(local_v,partition(gids)) |> consistent! |> wait
  local_v
end

## Parallel aggregation

using GridapEmbedded.Distributed: DistributedEmbeddedDiscretization

function _aggregate(
  stragegy,
  cutgeo::DistributedEmbeddedDiscretization,
  in_or_out::Integer=IN)

  map(local_views(cutgeo)) do cutgeo
    aggregate(stragegy,cutgeo,cutgeo.geo,in_or_out)
  end
end

# Driver

distribute = identity

np = (2,2)

ranks = distribute(LinearIndices((prod(np),)))

L = 1
p0 = Point(0.0,0.0)
pmin = p0-L/2
pmax = p0+L/2


R = 0.35
geo = disk(R,x0=p0)

n = 8
mesh_partition = (n,n)
bgmodel = CartesianDiscreteModel(ranks,np,pmin,pmax,mesh_partition)

cutgeo = cut(bgmodel,geo)

bgf_to_ioc = compute_bgfacet_to_inoutcut(bgmodel,geo)

Ω = Triangulation(cutgeo)

writevtk(Ω,"trian")

strategy = AggregateCutCellsByThreshold(0.8)
aggregates,aggregate_owner,aggregate_neig = _aggregate(strategy,cutgeo,geo,IN)

# TODO:
# - print aggregates
# - move to src
# - test aggregates with several geometries
# - reconstruct paths
# - add remote roots to model


# TODO: parallel aggregation (parallel aggfem article)
# 0. exchange measures and aggregates through gid
# 1. root in ghost layer
# 2. root in neighbors
# 3. root in neighbors of neighbors

end # module
