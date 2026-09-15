module NonConformingGridTopologies

using Gridap
using Gridap.Helpers

using Gridap.Geometry: GridTopology
using Gridap.Geometry: UnstructuredGridTopology
using Gridap.Geometry: get_grid_topology
using Gridap.Adaptivity: AdaptedDiscreteModel
import Gridap.Geometry: get_faces
using Gridap.Arrays: Table, length_to_ptrs!, lazy_append

using GridapDistributed

using GridapP4est: OctreeDistributedDiscreteModel
using GridapP4est: NonConformingGlue

using PartitionedArrays

struct NonConformingGridTopology{Dc,Dp,T,O} <: GridTopology{Dc,Dp}
  conforming_grid_topology::UnstructuredGridTopology{Dc,Dp,T,O}
  coarse_cell_to_hanging_faces::Vector{Table{Int32,Vector{Int32},Vector{Int32}}}
  hanging_faces_to_coarse_cells::Vector{Table{Int32,Vector{Int32},Vector{Int32}}}
  function NonConformingGridTopology(
    conforming_grid_topology::UnstructuredGridTopology{Dc,Dp,T,O},
    coarse_cell_to_hanging_faces::Vector{Table{Int32,Vector{Int32},Vector{Int32}}},
    hanging_faces_to_coarse_cells::Vector{Table{Int32,Vector{Int32},Vector{Int32}}}) where {Dc,Dp,T,O}
    @warn "NonConformingGridTopology is experimental. For now, it only provides the full grid 
    topology (conforming + non-conforming) when the dimfrom or dimto of the face to face maps 
    (get_faces) is Dc (dim of the cell)." maxlog=1
    new{Dc,Dp,T,O}(conforming_grid_topology,
                   coarse_cell_to_hanging_faces,
                   hanging_faces_to_coarse_cells)
  end
end

function NonConformingGridTopology(model::OctreeDistributedDiscreteModel)
  map(local_views(model),model.non_conforming_glue) do model,ncg
    conforming_grid_topology = get_grid_topology(model)
    coarse_cell_to_hanging_faces, hanging_faces_to_coarse_cells =
      generate_hanging_faces_to_coarse_cells_glue(model,ncg,
                                                 num_cells(model),
                                                 num_dims(model))
    NonConformingGridTopology(conforming_grid_topology,
                              coarse_cell_to_hanging_faces,
                              hanging_faces_to_coarse_cells)
  end
end

function generate_hanging_faces_to_coarse_cells_glue(
  model::DiscreteModel,ncg::NonConformingGlue,num_cells::Int,D::Int)

  num_cell_vertices, num_cell_edges, num_cell_faces =     
    _get_num_cell_vertices_edges_faces(model,D)
  topology = Gridap.Geometry.get_grid_topology(model)

  function _compute_hanging_face_to_cells(hanging_face_glue_tag)
    ocell, ocell_lface, _ = hanging_face_glue_tag
    if (ocell==-1)
      faces_to_cells = get_faces(topology,0,D)
      return similar(faces_to_cells[1], 0)
    end
    lface_dim = _face_dim(num_cell_vertices,
                          num_cell_edges,
                          num_cell_faces,
                          ocell_lface)
    ocell_lface_within_dim = _face_lid_within_dim(num_cell_vertices,
                                                  num_cell_edges,
                                                  num_cell_faces,
                                                  ocell_lface)
    cell_to_faces = get_faces(topology,D,lface_dim)
    owner_face = cell_to_faces[ocell][ocell_lface_within_dim]
    faces_to_cells = get_faces(topology,lface_dim,D)
    faces_to_cells[owner_face]
  end

  hanging_faces_to_num_cells, hanging_faces_to_cells_data,
  cell_to_hanging_face = 
    map(ncg.hanging_faces_glue) do glue
      data = map(_compute_hanging_face_to_cells, glue)
      lengths = map(length, data)
      total = sum(lengths)
      reduced = Vector{Int32}(undef, total)
      ids = Vector{Int32}(undef, total)

      pos = 1
      for i in eachindex(data)
        cells = data[i]
        m = length(cells)
        if m > 0
          copyto!(reduced, pos, cells, 1, m)
          ids[pos:pos+m-1] .= Int32(i)
          pos += m
        end
      end

      # optional safety check if you want to be strict:
      @assert pos - 1 == total
      lengths, reduced, ids
    end |> tuple_of_arrays

  hanging_faces_to_cells_ptrs = map(
      hanging_faces_to_num_cells,
      ncg.num_regular_faces,
      ncg.num_hanging_faces) do lengths, num_regular_faces, num_hanging_faces
    ptrs = zeros(Int32, num_regular_faces + num_hanging_faces + 1)
    ptrs[end-num_hanging_faces+1:end] .= lengths
    ptrs
  end
  map(length_to_ptrs!, hanging_faces_to_cells_ptrs)

  perm = map(sortperm, hanging_faces_to_cells_data)
  cell_to_hanging_faces_data = 
    map(cell_to_hanging_face,perm) do ids,p
      collect(Int32,lazy_map(Reindex(ids),p))
    end
  cell_to_hanging_faces_data = 
    map(Broadcasting(+),cell_to_hanging_faces_data,Int32.(ncg.num_regular_faces))

  cell_to_hanging_faces_ptrs = 
    map(_count_owned_hanging_faces,
        hanging_faces_to_cells_data,
        tfill(num_cells,Val{D}()))
  map(length_to_ptrs!,cell_to_hanging_faces_ptrs)
  
  coarse_cell_to_hanging_faces = map(cell_to_hanging_faces_data,
                                     cell_to_hanging_faces_ptrs) do d,p
    Table(d,p)
  end
  hanging_faces_to_coarse_cells = map(hanging_faces_to_cells_data,
                                      hanging_faces_to_cells_ptrs) do d,p
    Table(d,p)
  end

  coarse_cell_to_hanging_faces, hanging_faces_to_coarse_cells
end

function _get_num_cell_vertices_edges_faces(model,D::Int)
  cell_reffe = first(Gridap.Geometry.get_reffes(model))
  num_cell_vertices = Gridap.ReferenceFEs.num_faces(cell_reffe,0)
  num_cell_edges    = D==3 ? Gridap.ReferenceFEs.num_faces(cell_reffe,1) : 0
  num_cell_faces    = Gridap.ReferenceFEs.num_faces(cell_reffe,D-1)
  num_cell_vertices, num_cell_edges, num_cell_faces
end

function _face_dim(num_vertices,num_edges,num_faces,face_lid)
  if (face_lid <= num_vertices)
    return 0
  elseif (face_lid <= num_vertices + num_edges)
    return 1
  elseif (face_lid <= num_vertices + num_edges + num_faces)
    return (num_edges == 0 ? 1 : 2)
  end
end

function _face_lid_within_dim(num_vertices,num_edges,num_faces,face_lid)
  if (face_lid <= num_vertices)
    return face_lid
  elseif (face_lid <= num_vertices + num_edges)
    return face_lid - num_vertices
  elseif (face_lid <= num_vertices + num_edges + num_faces)
    return face_lid - num_vertices - num_edges
  end
end

function _count_owned_hanging_faces(x::AbstractArray{<:Integer},num_cells::Int)
  r = zeros(Int32,num_cells+1)
  length(x) == 0 && return r
  for xi in x
    r[xi+1] += 1
  end
  return r
end

export NonConformingGridTopology

# Implementation of abstract API

OrientationStyle(
  ::Type{<:NonConformingGridTopology{Dc,Dp,T,O}}) where 
  {Dc,Dp,T,O} = O()

RegularityStyle(
  ::Type{<:NonConformingGridTopology{Dc,Dp,T,O}}) where 
  {Dc,Dp,T,O} = Irregular()

get_vertex_coordinates(ncgt::NonConformingGridTopology) = 
  ncgt.conforming_grid_topology.vertex_coordinates

get_cell_type(ncgt::NonConformingGridTopology) = 
  ncgt.conforming_grid_topology.cell_type

get_polytopes(ncgt::NonConformingGridTopology) = 
  collect1d(ncgt.conforming_grid_topology.polytopes)

function get_faces(ncgt::NonConformingGridTopology{Dc,Dp,T,O}, 
                    dimfrom::Integer,dimto::Integer) where {Dc,Dp,T,O}
  get_faces(ncgt,Val{dimfrom}(),Val{dimto}())
end

function get_faces(ncgt::NonConformingGridTopology{Dc,Dp,T,O}, 
                    ::Val{Df},::Val{Dt}) where {Dc,Dp,T,O,Df,Dt}
  get_faces(ncgt.conforming_grid_topology,Df,Dt)
end

function get_faces(ncgt::NonConformingGridTopology{Dc,Dp,T,O}, 
                    ::Val{Dc},::Val{Dt}) where {Dc,Dp,T,O,Dt}
  r = get_faces(ncgt.conforming_grid_topology,Dc,Dt)
  h = ncgt.coarse_cell_to_hanging_faces[Dt+1]
  # Not sure which one of these is better
  # map(lazy_append,r,h)
  map(vcat,r,h)
end

function get_faces(ncgt::NonConformingGridTopology{Dc,Dp,T,O}, 
                    ::Val{Df},::Val{Dc}) where {Dc,Dp,T,O,Df}
  r = get_faces(ncgt.conforming_grid_topology,Df,Dc)
  h = ncgt.hanging_faces_to_coarse_cells[Df+1]
  # Not sure which one of these is better
  # map(lazy_append,r,h)
  map(vcat,r,h)
end

# Overload case dimfrom == dimto == Dc to avoid ambiguity
function get_faces(ncgt::NonConformingGridTopology{Dc,Dp,T,O}, 
                    ::Val{Dc},::Val{Dc}) where {Dc,Dp,T,O}
  get_faces(ncgt.conforming_grid_topology,Dc,Dc)
end

end