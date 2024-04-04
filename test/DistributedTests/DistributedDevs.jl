module DistributedDevs

using Gridap
using GridapEmbedded
using GridapDistributed
using PartitionedArrays
using Test

using GridapEmbedded.CSG
using GridapEmbedded.Interfaces: CUT
using Gridap.Geometry
using Gridap.ReferenceFEs

using GridapEmbedded.AgFEM: _touch_aggregated_cells!

function _aggregate_one_step!(
  c1,c2,own_cells,
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

        cell_to_neig[cell] = neigh_cell
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

function _aggregate(
  stragegy,
  cutgeo::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out::Integer=IN)

  bgf_to_ioc = compute_bgfacet_to_inoutcut(bgmodel,geo)
  _aggregate(stragegy,cutgeo,geo,in_or_out,bgf_to_ioc)
end

function _aggregate(
  stragegy,
  cutgeo::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out::Integer,
  bgf_to_ioc::AbstractArray{<:AbstractVector})

  map(local_views(cutgeo),bgf_to_ioc) do cutgeo,bgf_to_ioc
    #input cell meas or IN cells
    aggregate(stragegy,cutgeo,geo,in_or_out,bgf_to_ioc)
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





threshold = 0.8
loc = IN
facet_to_inoutcut = compute_bgfacet_to_inoutcut(model,geo)


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

ocell_to_threshold = map(cell_to_unit_cut_meas) do c_to_m
  n_cells = length(c_to_m)
  c_to_t = falses(n_cells)
  c = array_cache(c_to_m)
  for cell in 1:n_cells
    if getindex!(c,c_to_m,cell) ≥ threshold
      c_to_t[cell] = true
    end
  end
  c_to_t
end

cell_to_root_centroid = map(cell_to_coords) do cell_to_coords
  map(i->sum(i)/length(i),cell_to_coords)
end
PVector(cell_to_root_centroid,partition(gids)) |> consistent! |> wait


gids = get_cell_gids(model)

cell_to_threshold = map(ocell_to_threshold,local_to_own(gids)) do o_to_t,l_to_o
  T = eltype(o_to_t)
  map(o-> iszero(o) ? zero(T) : o_to_t[o],l_to_o)
end
PVector(cell_to_threshold,partition(gids)) |> consistent! |> wait
cell_to_threshold

n_cells = map(length,cell_to_threshold)
cell_to_cellin = map(n->zeros(Int32,n),n_cells)
cell_to_touched = map(falses,n_cells)
cells_threshold = map(findall,cell_to_threshold)
cell_to_cellin = map(cells_threshold,cell_to_cellin) do cells,c_to_ci
  c_to_ci[cells] = cells # TODO: gid here
  c_to_ci
end
cell_to_neig = map(copy,cell_to_cellin)
cell_to_touched = cell_to_threshold


c1 = map(array_cache,cell_to_faces)
c2 = map(array_cache,face_to_cells)

max_iters = 1
all_aggregated = false

# for iter in 1:max_iters
own_cells = own_to_local(gids)[1]

cell_to_root_part = map(collect,local_to_owner(gids))
# communicate

# end # for


# init variables


# for iter in 1:max_iters

all_aggregated = map(c1,c2,own_to_local(gids),
  cell_to_inoutcut,
  cell_to_touched,
  cell_to_neig,
  cell_to_cellin,
  cell_to_root_centroid,
  cell_to_root_part,
  cell_to_faces,
  face_to_cells,
  facet_to_inoutcut) do c1,c2,own_cells,
      cell_to_inoutcut,
      cell_to_touched,
      cell_to_neig,
      cell_to_cellin,
      cell_to_root_centroid,
      cell_to_root_part,
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
    cell_to_faces,
    face_to_cells,
    facet_to_inoutcut,
    loc)

end


# consistent! (communicate ghost)
(cell_to_touched,
cell_to_neig, # neig_lid (if lid not consistent, ignore ghost)
cell_to_cellin, # root_gid
cell_to_root_centroid,
cell_to_root_part)

# reduce all_aggregated and break if all(all_aggregated)


# return cell_to_rood_gid, cell_to_root_part, cell_to_neig

# BUG: aggregates do not have information of ghost cut cells
# TODO: exchange measures and aggregates
# _aggregate(strategy,cutgeo)


# TODO: parallel aggregation (parallel aggfem article)
# 0. exchange measures and aggregates through gid
# 1. root in ghost layer
# 2. root in neighbors
# 3. root in neighbors of neighbors

end # module
