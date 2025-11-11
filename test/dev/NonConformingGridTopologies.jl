module NonConformingGridTopologies

  using Gridap
  using Gridap.Helpers

  using Gridap.Geometry: GridTopology
  using Gridap.Geometry: UnstructuredGridTopology
  using Gridap.Geometry: get_grid_topology
  import Gridap.Geometry: get_faces
  using Gridap.Arrays: Table, length_to_ptrs!, lazy_append

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
      map(Broadcasting(zeros),tfill(Int32,Val{D}()),
          ncg.num_regular_faces.+ncg.num_hanging_faces.+1)
    map(_ones_on_hanging_faces!,hanging_faces_to_cell_ptrs,
                                ncg.num_hanging_faces)
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

    coarse_cell_to_hanging_faces, hanging_faces_to_coarse_cell
  end

  function _ones_on_hanging_faces!(
    ptrs::Vector{Int32},num_hanging_faces::Int)
    ptrs[end-num_hanging_faces+1:end] .= 1
    return nothing
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
    h = ncgt.hanging_faces_to_coarse_cell[Df+1]
    # Not sure which one of these is better
    # map(lazy_append,r,h)
    map(vcat,r,h)
  end

  # Overload case dimfrom == dimto == Dc to avoid ambiguity
  function get_faces(ncgt::NonConformingGridTopology{Dc,Dp,T,O}, 
                     ::Val{Dc},::Val{Dc}) where {Dc,Dp,T,O}
    get_faces(ncgt.conforming_grid_topology,Dc,Dc)
  end

  # What other procs do I need? Review CA (e.g., facet_to_inoutcut)

end # module