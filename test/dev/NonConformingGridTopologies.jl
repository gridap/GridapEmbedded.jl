module NonConformingGridTopologies

  using Gridap
  using Gridap.Helpers

  import Gridap.Geometry: get_faces
  using Gridap.Geometry: GridTopology
  using Gridap.Geometry: UnstructuredGridTopology
  using Gridap.Geometry: get_grid_topology
  using Gridap.Arrays: Table, length_to_ptrs!

  using GridapDistributed

  using GridapP4est: OctreeDistributedDiscreteModel
  using GridapP4est: NonConformingGlue

  struct NonConformingGridTopology{Dc,Dp,T,O} <: GridTopology{Dc,Dp}
    conforming_grid_topology::UnstructuredGridTopology{Dc,Dp,T,O}
    coarse_cell_to_hanging_faces::Vector{Table{Int32,Vector{Int32},Vector{Int32}}}
    hanging_faces_to_coarse_cell::Vector{Table{Int32,Vector{Int32},Vector{Int32}}}
  end

  function NonConformingGridTopology(model::OctreeDistributedDiscreteModel)
    map(local_views(model),model.non_conforming_glue) do model,ncg
      conforming_grid_topology = get_grid_topology(model)
      coarse_cell_to_hanging_faces, hanging_faces_to_coarse_cell =
        generate_hanging_faces_to_coarse_cell_glue(ncg,
                                                   num_cells(model),
                                                   num_dims(model))
      NonConformingGridTopology(conforming_grid_topology,
                                coarse_cell_to_hanging_faces,
                                hanging_faces_to_coarse_cell)
    end
  end

  function generate_hanging_faces_to_coarse_cell_glue(
    ncg::NonConformingGlue,num_cells::Int,D::Int)

    hanging_faces_to_cell_data = 
      map(Broadcasting(first),ncg.hanging_faces_glue)

    hanging_faces_to_cell_ptrs = 
      map(Broadcasting(ones),tfill(Int32,Val{D}()),ncg.num_hanging_faces.+1)
    map(length_to_ptrs!,hanging_faces_to_cell_ptrs)
    
    cell_to_hanging_faces_data = 
      map(sortperm,hanging_faces_to_cell_data)
    cell_to_hanging_faces_data =
      map(Broadcasting(+),cell_to_hanging_faces_data,ncg.num_regular_faces)
    
    cell_to_hanging_faces_ptrs = 
      map(_count_owned_hanging_faces,
          hanging_faces_to_cell_data,
          tfill(num_cells,Val{D}()))
    map(length_to_ptrs!,cell_to_hanging_faces_ptrs)
    
    # Convert data from Int64 to Int32
    cell_to_hanging_faces_data =
      map(Broadcasting(Int32),cell_to_hanging_faces_data)
    hanging_faces_to_cell_data =
      map(Broadcasting(Int32),hanging_faces_to_cell_data)

    coarse_cell_to_hanging_faces = map(cell_to_hanging_faces_data,
                                       cell_to_hanging_faces_ptrs) do d,p
      Table(d,p)
    end
    hanging_faces_to_coarse_cell = map(hanging_faces_to_cell_data,
                                       hanging_faces_to_cell_ptrs) do d,p
      Table(d,p)
    end

    hanging_faces_to_coarse_cell, coarse_cell_to_hanging_faces
  end

  function _count_owned_hanging_faces(x::AbstractArray{<:Integer},num_cells::Int)
    r = zeros(Int32,num_cells+1)
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
                     ::Val{Dc},dimto::Integer) where {Dc,Dp,T,O}
    @show ncgt.coarse_cell_to_hanging_faces[dimto+1]
  end

  function get_faces(ncgt::NonConformingGridTopology{Dc,Dp,T,O}, 
                     dimfrom::Integer,::Val{Dc}) where {Dc,Dp,T,O}
    @show ncgt.hanging_faces_to_coarse_cell[dimfrom+1]
  end

  # What other procs do I need? Review CA (e.g., facet_to_inoutcut)

end # module