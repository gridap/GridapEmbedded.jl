module DistributedAggregationP4estMeshes

  using Gridap
  using GridapEmbedded
  using GridapDistributed
  using PartitionedArrays
  using MPI

  using P4est_wrapper
  using GridapP4est

  using Gridap.Arrays: Table
  using Gridap.Geometry: get_grid_topology
  using Gridap.Geometry: get_faces

  using FillArrays

  function run(distribute)

    ranks = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))
    coarse_model = CartesianDiscreteModel((-1,1,-1,1,-1,1),(1,1,1))
    num_uniform_refinements = 1
    num_ghost_layers = 2
    D = num_dims(coarse_model)

    dmodel = OctreeDistributedDiscreteModel(ranks,
                                            coarse_model,
                                            num_uniform_refinements;
                                            num_ghost_layers=num_ghost_layers)
    
    fmodel_refine_coarsen_flags = 
      map(ranks,partition(get_cell_gids(dmodel.dmodel))) do rank,indices
        flags = zeros(Int,length(indices))
        flags .= nothing_flag        
        flags[1] = refine_flag
        flags
    end
    fmodel,_ = Gridap.Adaptivity.adapt(dmodel,fmodel_refine_coarsen_flags);

    # map(local_views(fmodel)) do m
    #   gt = get_grid_topology(m)
    #   @show get_faces(gt,1,2)
    #   @show get_faces(gt,2,1)
    # end

    # map(fmodel.non_conforming_glue) do ncg
    #   @show ncg.num_regular_faces
    #   @show ncg.num_hanging_faces
    #   @show ncg.hanging_faces_glue
    #   @show ncg.hanging_faces_to_cell
    #   @show ncg.hanging_faces_to_lface
    #   @show ncg.owner_faces_pindex
    #   @show ncg.owner_faces_lids
    # end

    cell_to_hanging_faces = 
      map(local_views(fmodel),fmodel.non_conforming_glue) do m,ncg
      
      hanging_face_to_owner_cell = 
        map(Broadcasting(first),ncg.hanging_faces_glue)
      cell_to_hanging_faces_data = 
        map(sortperm,hanging_face_to_owner_cell)
      @time hanging_face_to_owner_cell = 
        map(sort,hanging_face_to_owner_cell)

      # # Only for facets
      # hanging_face_to_owner_cell = first.(ncg.hanging_faces_glue[D])
      # cell_to_hanging_faces_data = sortperm(hanging_face_to_owner_cell)

      # cell_to_hanging_faces_ptr = zeros(Int32,num_cells(m)+1)
      # owner_cells_ptr = unique(hanging_face_to_owner_cell) .+ 1
      # # Alternatively, for owner_cells_ptr evaluate every 2^(D-1) entries
      # # But with unique I can do it for every dimension
      # cell_to_hanging_faces_ptr[owner_cells_ptr] .= 2^(D-1)
      # length_to_ptrs!(cell_to_hanging_faces_ptr)

      # cell_to_hanging_faces_ptr = [zeros(Int32,num_cells(m)+1) for i in 1:D] 

      # # Non-sorted 3D: 0.028385 seconds (65.36 k allocations: 3.451 MiB, 99.81% compilation time)
      # # Sorted 3D: 0.027257 seconds (65.36 k allocations: 3.451 MiB, 99.73% compilation time)
      # @time map(cell_to_hanging_faces_ptr,hanging_face_to_owner_cell) do cthf_ptr,hftoc
      #   for oc in hftoc
      #     cthf_ptr[oc+1] += 1
      #   end
      # end
      @time owner_cells = map(unique,hanging_face_to_owner_cell)
      # # Sorted 3D: 0.039553 seconds (113.76 k allocations: 6.096 MiB, 99.85% compilation time)
      # @time cell_to_hanging_faces_ptr = map(hanging_face_to_owner_cell,owner_cells) do hftoc,oc
      #   map(x->searchsortedfirst(hftoc,x),oc)
      # end
      # # Sorted 3D: 0.000010 seconds (22 allocations: 1.766 KiB)
      @time owner_cell_to_hanging_faces_ptr = 
        map(indexin,owner_cells,hanging_face_to_owner_cell)

      cell_to_hanging_faces_data =
        map(Broadcasting(+),cell_to_hanging_faces_data,ncg.num_regular_faces)
      # cell_to_hanging_faces_data = 
      #   cell_to_hanging_faces_data .+ ncg.num_regular_faces[D]

      # @show cell_to_hanging_faces_ptr
      @show hanging_face_to_owner_cell
      @show cell_to_hanging_faces_data
      @show owner_cell_to_hanging_faces_ptr
      @show owner_cells

      # Table(cell_to_hanging_faces_ptr,cell_to_hanging_faces_data)
      cell_to_hanging_faces_data
    end

    # @show typeof(cell_to_hanging_faces)

  end

end