module DistributedAgFEM4P4estMeshes

using Gridap
using Gridap.Helpers
using GridapEmbedded
using GridapEmbedded.Interfaces: CUT
using GridapEmbedded.LevelSetCutters: Leaf
using GridapEmbedded.Distributed: _local_aggregates
using GridapDistributed
using PartitionedArrays
using MPI

using Gridap.Geometry: get_faces
using Gridap.FESpaces
using Gridap.Arrays

using P4est_wrapper
using GridapP4est
using GridapP4est: generate_local_fe_spaces_and_constraints

using GridapDistributed: DistributedCellField, DistributedCellDof

using DrWatson
using Plots

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

function find_optimal_roots!(lcell_to_root,lcell_to_value,lcell_to_owner,cell_indices)
  # Bring all root candidates to the owner of the cut cell
  roots_cache = PartitionedArrays.p_vector_cache(lcell_to_root,cell_indices)
  t1 = exchange_impl!(lcell_to_root,roots_cache)
  values_cache = PartitionedArrays.p_vector_cache(lcell_to_value,cell_indices)
  t2 = exchange_impl!(lcell_to_value,values_cache)
  owners_cache = PartitionedArrays.p_vector_cache(lcell_to_owner,cell_indices)
  t3 = exchange_impl!(lcell_to_owner,owners_cache)
  wait(t1)
  wait(t2)
  wait(t3)

  # Select the optimal root for each local cut cell
  map(
    lcell_to_root,lcell_to_value,lcell_to_owner,roots_cache,values_cache,owners_cache
  ) do lcell_to_root, lcell_to_value, lcell_to_owner, roots_cache, values_cache, owners_cache
    lids_rcv = roots_cache.local_indices_rcv
    values_rcv = values_cache.buffer_rcv
    roots_rcv = roots_cache.buffer_rcv
    owners_rcv = owners_cache.buffer_rcv
    for k in eachindex(roots_cache.neighbors_rcv)
      for (lcell,root,value,owner) in zip(lids_rcv[k], roots_rcv[k], values_rcv[k], owners_rcv[k])
        if lcell_to_value[lcell] > value # Take the minimum
          lcell_to_root[lcell] = root # Global ID
          lcell_to_value[lcell] = value
          lcell_to_owner[lcell] = owner
        end
      end
    end
  end

  # Scatter the optimal roots and values
  # Technically, we do not need to scatter the values, it can be removed after 
  # we are done debugging
  t1 = consistent!(PVector(lcell_to_owner,cell_indices,owners_cache))
  t2 = consistent!(PVector(lcell_to_root,cell_indices,roots_cache))
  t3 = consistent!(PVector(lcell_to_value,cell_indices,values_cache))
  
  wait(t1)

  agg_cell_indices = map(cell_indices,lcell_to_owner) do cell_indices, lcell_to_owner
    LocalIndices(global_length(cell_indices),part_id(cell_indices),local_to_global(cell_indices),lcell_to_owner)
  end

  wait(t2)
  wait(t3)

  return lcell_to_root, lcell_to_owner, lcell_to_value, agg_cell_indices
end

# This respects the input triangulation, without adding back all the ghosts
function MyFESpace(
  trian::GridapDistributed.DistributedTriangulation,reffe,spaces,
  n_fdofs,n_ddofs,gids,
  mDOF_to_dof,sDOF_to_dof,sDOF_to_mdofs,sDOF_to_coeffs;
  split_own_and_ghost=false,constraint=nothing,kwargs...
)
  myfespaces = map(mDOF_to_dof,sDOF_to_dof,sDOF_to_mdofs,sDOF_to_coeffs,spaces,n_fdofs,n_ddofs) do mDOF_to_dof,sDOF_to_dof,sDOF_to_mdofs,sDOF_to_coeffs,spaces,n_fdofs,n_ddofs
    FESpaceWithLinearConstraints(spaces,
                                 mDOF_to_dof,
                                 sDOF_to_dof,
                                 sDOF_to_mdofs,
                                 sDOF_to_coeffs,
                                 n_fdofs = n_fdofs,
                                 n_dofs = n_fdofs + n_ddofs)
  end
  vector_type = GridapDistributed._find_vector_type(myfespaces,gids;split_own_and_ghost=split_own_and_ghost)
  space = GridapDistributed.DistributedSingleFieldFESpace(myfespaces,gids,trian,vector_type)
  return GridapDistributed._add_distributed_constraint(space,reffe,constraint)
end

function Gridap.Adaptivity.get_model(model::GridapP4est.OctreeDistributedDiscreteModel)
  dmodel = GridapDistributed.GenericDistributedDiscreteModel(
    map(Gridap.Adaptivity.get_model,local_views(model)),
    get_cell_gids(model)
  )
  return OctreeDistributedDiscreteModel(
    num_cell_dims(model), num_point_dims(model),
    model.parts,dmodel,model.non_conforming_glue,
    model.coarse_model,
    model.ptr_pXest_connectivity,
    model.ptr_pXest,
    model.pXest_type,
    model.pXest_refinement_rule_type,
    model.connectivity_ref; num_ghost_layers = model.num_ghost_layers
  )
end

include("NonConformingGridTopologies.jl")
using .NonConformingGridTopologies

abstract type UnfittedGeometry end

function generate_geometry(::UnfittedGeometry)
  @notimplemented "`generate_geometry` must be implemented for this UnfittedGeometry subtype"
end

struct Square <: UnfittedGeometry end

struct SquareWithCoarseCuts <: UnfittedGeometry end

function generate_geometry(::Union{Square,SquareWithCoarseCuts})
  quadrilateral(;x0=Point(-0.88,-0.88),
                 d1=VectorValue(1.75,0.0),
                 d2=VectorValue(0.0,1.75))
  # quadrilateral(;x0=Point(-0.750001,-0.750001),
  #                d1=VectorValue(1.500002,0.0),
  #                d2=VectorValue(0.0,1.500002))
end

struct SquareWithCircularHole <: UnfittedGeometry end

function generate_geometry(::SquareWithCircularHole)
  geo1 = quadrilateral(;x0=Point(-0.9,-0.9),
                        d1=VectorValue(1.8,0.0),
                        d2=VectorValue(0.0,1.8))
  geo2 = ! disk(0.4)
  geo = intersect(geo1,geo2)
  return geo
end

struct Flower <: UnfittedGeometry
  x₀
  R₀::Float64
  m::Float64
  ω::Float64
end

Flower(;x₀=Point(0.0,0.0),R₀=0.6,m=0.6,ω=10.0) = Flower(x₀,R₀,m,ω)

function generate_geometry(geometry::Flower)
  name="flower"
  function flowerfun(x)
    _flower(x,geometry.x₀,geometry.R₀,geometry.m,geometry.ω)
  end
  tree = Leaf((flowerfun,name,nothing))
  geo = AnalyticalGeometry(tree)
  return geo
end

@inline function _flower(x::Point,x₀,R₀,m,ω)
  w = x - x₀
  t = angle(w[1]+w[2]*im)
  w⋅w - (R₀*(1.0+m*sin(ω*t)))^2
end

struct SinusoidalBand <: UnfittedGeometry
  y₀::Float64
  thickness::Float64
  lower_amplitude::Float64
  lower_period::Float64
  upper_amplitude::Float64
  upper_period::Float64
  lower_phase::Float64
  upper_phase::Float64
end

SinusoidalBand(;y₀=0.0,
                thickness=0.6,
                lower_amplitude=0.08,
                lower_period=1.1,
                upper_amplitude=0.12,
                upper_period=1.7,
                lower_phase=0.0,
                upper_phase=pi/3) =
  SinusoidalBand(y₀,
                 thickness,
                 lower_amplitude,
                 lower_period,
                 upper_amplitude,
                 upper_period,
                 lower_phase,
                 upper_phase)

function generate_geometry(geometry::SinusoidalBand)
  name="sinusoidal_band"
  function bandfun(x)
    _sinusoidal_band(x,geometry)
  end
  tree = Leaf((bandfun,name,nothing))
  geo = AnalyticalGeometry(tree)
  return geo
end

@inline function _sinusoidal_band(x::Point,geometry::SinusoidalBand)
  ξ = x[1]
  y = x[2]
  lower = geometry.y₀ - geometry.thickness/2 +
          geometry.lower_amplitude*sin(2*pi*ξ/geometry.lower_period + geometry.lower_phase)
  upper = geometry.y₀ + geometry.thickness/2 +
          geometry.upper_amplitude*sin(2*pi*ξ/geometry.upper_period + geometry.upper_phase)
  return max(lower - y, y - upper)
end

function generate_unfitted_model(ranks,
                                 num_parts,
                                 geometry::UnfittedGeometry;
                                 initial_uniform_refs=4,
                                 post_uniform_refinements=0,
                                 num_ghost_layers::Integer=1)
  coarse_model = CartesianDiscreteModel((-1,1,-1,1),(1,1))
  dmodel = OctreeDistributedDiscreteModel(
    ranks,coarse_model,initial_uniform_refs;num_ghost_layers=num_ghost_layers
  )
  geo = generate_geometry(geometry)
  cutgeo = cut(dmodel,geo)
  cell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geo)
  fmodel_refine_coarsen_flags = 
    map(ranks,
        partition(get_cell_gids(dmodel.dmodel)),
        cell_to_inoutcut) do rank,indices,cell_to_inoutcut
      flags = zeros(Int,length(indices))
      flags .= nothing_flag
      toref = findall(c->c!=CUT,cell_to_inoutcut)
      flags[toref] .= coarsen_flag
      flags
  end
  fmodel,_ = Gridap.Adaptivity.adapt(dmodel,fmodel_refine_coarsen_flags);
  for i in 1:initial_uniform_refs-1
    cutgeo = cut(fmodel,geo)
    cell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geo)
    fmodel_refine_coarsen_flags = 
      map(partition(get_cell_gids(fmodel)),
          cell_to_inoutcut) do indices,cell_to_inoutcut
        flags = zeros(Int,length(indices))
        flags .= nothing_flag
        toref = findall(c->c!=CUT,cell_to_inoutcut)
        flags[toref] .= coarsen_flag
        flags
    end
    fmodel,_ = Gridap.Adaptivity.adapt(fmodel,fmodel_refine_coarsen_flags);
  end
  if ( geometry isa SquareWithCoarseCuts ) && ( initial_uniform_refs == 4 )
    cutgeo = cut(fmodel,geo)
    fmodel_refine_coarsen_flags = 
      map(partition(get_cell_gids(fmodel))) do indices
        flags = zeros(Int,length(indices))
        flags .= nothing_flag
        flags[11:16] .= coarsen_flag
        flags[length(flags)-13:length(flags)-10] .= coarsen_flag
        flags
    end
    fmodel,_ = Gridap.Adaptivity.adapt(fmodel,fmodel_refine_coarsen_flags);
  end
  for i in 1:post_uniform_refinements
    cutgeo = cut(fmodel,geo)
    fmodel_refine_coarsen_flags = 
      map(partition(get_cell_gids(fmodel))) do indices
        flags = zeros(Int,length(indices))
        flags .= refine_flag
        flags
    end
    fmodel,_ = Gridap.Adaptivity.adapt(fmodel,fmodel_refine_coarsen_flags);
  end
  fmodel, _ = GridapDistributed.redistribute(fmodel)
  cutgeo = cut(fmodel,geo)
  return fmodel, cutgeo, geo
end

function generate_amr_constraints(_dtrian, spaces, reffe)
  model = get_background_model(_dtrian)
  trian = GridapDistributed.add_ghost_cells(_dtrian)
  cell_gids = GridapDistributed.generate_cell_gids(trian)
  models, non_conforming_glue = GridapP4est._generate_active_models_and_non_conforming_glue(
    model.pXest_type, model.pXest_refinement_rule_type, trian, cell_gids, model.non_conforming_glue
  )
  sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = GridapP4est.generate_constraints(
    model.pXest_refinement_rule_type, reffe, models, non_conforming_glue, spaces
  )
  # [TODO] amr dofs are Int64, agg dofs are Int32, we must unify this
  sDOF_to_dofs = map(sDOF_to_dofs) do sDOF_to_dofs
    Table(Int32.(sDOF_to_dofs.data), sDOF_to_dofs.ptrs)
  end
  return sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs
end

function generate_aggregates(fmodel, cutgeo, geo)
  cell_gids = get_cell_gids(fmodel)
  cell_indices = partition(cell_gids)

  ncgt = NonConformingGridTopology(fmodel)
  strategy = AggregateCutCellsByThreshold(1.0)
  lcell_to_lroot, lcell_to_root, lcell_to_value = map(
    local_views(cutgeo),cell_indices,ncgt
  ) do cutgeo,cell_indices,ncgt
    lid_to_gid = local_to_global(cell_indices)
    aggregate(strategy,cutgeo,geo,lid_to_gid,IN,grid_topology=ncgt)
  end |> tuple_of_arrays
  
  lcell_to_owner = map(copy∘local_to_owner,cell_indices)
  lcell_to_owner = map(lcell_to_owner,lcell_to_lroot) do lcell_to_owner,lcell_to_lroot
    for i in eachindex(lcell_to_owner)
      if !iszero(lcell_to_lroot[i])
        lcell_to_owner[i] = lcell_to_owner[lcell_to_lroot[i]]
      end
    end
    lcell_to_owner
  end

  lcell_to_root, lcell_to_owner, lcell_to_value, agg_cell_indices =
    find_optimal_roots!(lcell_to_root,lcell_to_value,lcell_to_owner,cell_indices);

  return lcell_to_root, agg_cell_indices
end

function generate_agfem_constraints(trian, spaces, bgcell_to_groot)
  shfns = DistributedCellField(map(get_fe_basis, spaces), trian)
  dofs = DistributedCellDof(map(get_fe_dof_basis, spaces), trian)
  model = get_background_model(trian)
  cell_gids = get_cell_gids(model)
  bgcell_to_gcell = local_to_global(cell_gids)
  # (TMP) This gives wrong output when root is not owned, but does 
  # not affect below, because the bgcell_to_bgroot is only queried
  # on owned cells
  bgcell_to_bgroot = _local_aggregates(bgcell_to_groot,cell_gids)
  sDOF_gids, new_n_fdofs, new_n_ddofs, mfdof_gids, mddof_gids, mDOF_to_dof, sDOF_to_dof, sDOF_to_DOFs, sDOF_to_coeffs = 
  GridapEmbedded.Distributed.generate_aggregated_space_constraints(
    trian, spaces, bgcell_to_bgroot, shfns, dofs, bgcell_to_gcell
  )
  return sDOF_gids, new_n_fdofs, new_n_ddofs, mfdof_gids, mddof_gids, mDOF_to_dof, sDOF_to_dof, sDOF_to_DOFs, sDOF_to_coeffs
end

function generate_merged_constraints(trian, 
                                     spaces, 
                                     bgcell_to_groot, 
                                     amr_sDOF_to_dof, 
                                     amr_sDOF_to_dofs, 
                                     amr_sDOF_to_coeffs)
  shfns = DistributedCellField(map(get_fe_basis, spaces), trian)
  dofs = DistributedCellDof(map(get_fe_dof_basis, spaces), trian)
  model = get_background_model(trian)
  cell_gids = get_cell_gids(model)
  bgcell_to_gcell = local_to_global(cell_gids)
  map(bgcell_to_groot,partition(cell_gids)) do bgcell_to_groot,cell_ids
    own_gids = local_to_global(cell_ids)[own_to_local(cell_ids)]
    owned_gids = Set(own_gids)
    for bgcell in own_to_local(cell_ids)
      groot = bgcell_to_groot[bgcell]
      groot == 0 && continue
      @assert groot in owned_gids "The root of an owned cell must also be owned"
    end
  end
  # (TMP) This gives wrong output when root is not owned, but does 
  # not affect below, because the bgcell_to_bgroot is only queried
  # on owned cells
  bgcell_to_bgroot = _local_aggregates(bgcell_to_groot,cell_gids)
  sDOF_gids, new_n_fdofs, new_n_ddofs, mfdof_gids, mddof_gids, mDOF_to_dof, sDOF_to_dof, sDOF_to_DOFs, sDOF_to_coeffs = 
    GridapEmbedded.Distributed.generate_merged_space_constraints(
      trian, 
      spaces, 
      bgcell_to_bgroot, 
      shfns, 
      dofs, 
      bgcell_to_gcell, 
      amr_sDOF_to_dof, 
      amr_sDOF_to_dofs, 
      amr_sDOF_to_coeffs
    )
  return sDOF_gids, new_n_fdofs, new_n_ddofs, mfdof_gids, mddof_gids, mDOF_to_dof, sDOF_to_dof, sDOF_to_DOFs, sDOF_to_coeffs
end

# ---------------------------------------------------------------------------
# Aggregate-conforming model with a single ghost layer
#
# Adapted from change_ownership_snippet.jl (change_ownership + trim_ghost_layers).
# generate_aggregates already gives the new owners and the cell indices
# (agg_cell_indices), so the snippet's step 1 is not needed here:
#  a) Wrap the local models and non-conforming glue of the input model,
#     unchanged, with agg_cell_indices. No geometry is rebuilt.
#  b) Drop the ghost cells beyond one layer of the new owned cells, using
#     Gridap's DiscreteModelPortion plus a restriction of the non-conforming glue.
# A part can only own a cell that it already holds, i.e., a cell in its
# input ghost layers. This is checked.
# The result is not backed by a p4est forest (ptr_pXest === nothing), since p4est
# cannot represent this partition. FE spaces, triangulations and integration
# work on it, but adapt/refine/coarsen/redistribute do not.
# ---------------------------------------------------------------------------

"""
    aggregate_conforming_model(model,agg_cell_indices) -> agg_model, new_to_old_cell

Returns an `OctreeDistributedDiscreteModel` with the same mesh as `model`, with
one ghost layer, in which the cells are owned as in `agg_cell_indices`. Global
cell ids are preserved. Also returns, for each part, the local ids in `model` of
the local cells of `agg_model`, to move cell-wise data (e.g., the aggregates) to
`agg_model`.
"""
function aggregate_conforming_model(model::OctreeDistributedDiscreteModel,agg_cell_indices)
  # Every cell must be owned by some part. It is not if its new owner does not
  # hold it as a local cell.
  num_owned_cells = reduction(+,map(own_length,agg_cell_indices);init=0,destination=:all)
  @assert PartitionedArrays.getany(num_owned_cells) == num_cells(model) """
    Some cells are given to a part that does not hold them.
    Increase the number of ghost layers of the input model."""

  # Number of ghost layers that the local cells of `model` guarantee around the
  # new owned cells. Each part needs at least one, i.e., it must hold every
  # cell touching a cell it owns.
  num_layers = map(local_views(model),partition(get_cell_gids(model)),agg_cell_indices) do lmodel,ids,agg_ids
    old_layer = cell_layers(lmodel,own_to_local(ids),model.num_ghost_layers)
    new_owned_layer = old_layer[own_to_local(agg_ids)]
    any(==(-1),new_owned_layer) && return 0
    model.num_ghost_layers - maximum(new_owned_layer;init=0)
  end
  agg_num_ghost_layers = PartitionedArrays.getany(reduction(min,num_layers;init=typemax(Int),destination=:all))
  @assert agg_num_ghost_layers >= 1 """
    Some part would own a cell in the outermost ghost layer of the input model.
    Increase the number of ghost layers of the input model."""

  # a) Same local cells, aggregate-conforming owners
  wrapped_model = octree_model_with_local_models(model,
                                                 local_views(model),
                                                 model.non_conforming_glue,
                                                 PRange(agg_cell_indices),
                                                 agg_num_ghost_layers)
  # b) Keep only one ghost layer
  return trim_ghost_layers(wrapped_model,1)
end

"""
    trim_ghost_layers(model,num_ghost_layers) -> new_model, new_to_old_cell

Restriction of `model` to its owned cells plus `num_ghost_layers` layers of ghost
cells. Uses Gridap's DiscreteModelPortion for the grid, topology and face
labeling, and restricts the non-conforming glue. Also returns, for each part,
the local ids in `model` of the local cells of `new_model`.
"""
function trim_ghost_layers(model::OctreeDistributedDiscreteModel,num_ghost_layers::Integer)
  @assert 1 <= num_ghost_layers <= model.num_ghost_layers
  ids = partition(get_cell_gids(model))

  # Kept cells (owned cells first), in the local numbering of model
  new_to_old_cell = map(local_views(model),ids) do lmodel,ids
    layer = cell_layers(lmodel,own_to_local(ids),num_ghost_layers)
    l2g   = local_to_global(ids)
    own   = sort(collect(own_to_local(ids)),by=c->l2g[c])
    gho   = sort(findall(>(0),layer),by=c->l2g[c])
    vcat(own,gho)
  end
  new_ids = map(ids,new_to_old_cell) do ids,n2o
    LocalIndices(global_length(ids),part_id(ids),
                 collect(Int,local_to_global(ids)[n2o]),
                 collect(Int32,local_to_owner(ids)[n2o]))
  end
  # The ghost relation is symmetric (cell distance), as in setup_cell_prange
  assembly_neighbors(new_ids;symmetric=true)

  models, non_conforming_glue = map(local_views(model),model.non_conforming_glue,new_to_old_cell) do lmodel,ncglue,n2o
    portion = Gridap.Geometry.DiscreteModelPortion(lmodel,n2o)
    portion.model, restrict_non_conforming_glue(model.pXest_refinement_rule_type,ncglue,portion)
  end |> tuple_of_arrays
  new_model = octree_model_with_local_models(model,models,non_conforming_glue,PRange(new_ids),Int(num_ghost_layers))
  return new_model, new_to_old_cell
end

# OctreeDistributedDiscreteModel with the same coarse model, connectivity and
# refinement rule as `model`, but new local models and cell gids. It is not
# backed by a p4est forest, as p4est cannot represent arbitrary partitions.
# (Passing model.ptr_pXest instead would make both models destroy the same
# forest when they are garbage collected.)
function octree_model_with_local_models(model::OctreeDistributedDiscreteModel{Dc,Dp},
                                        models,non_conforming_glue,cell_gids,
                                        num_ghost_layers) where {Dc,Dp}
  dmodel = GridapDistributed.DistributedDiscreteModel(models,cell_gids)
  OctreeDistributedDiscreteModel(Dc,Dp,model.parts,dmodel,non_conforming_glue,
                                 model.coarse_model,
                                 model.ptr_pXest_connectivity,
                                 nothing,
                                 model.pXest_type,
                                 model.pXest_refinement_rule_type,
                                 model.connectivity_ref;
                                 num_ghost_layers=num_ghost_layers)
end

# Distance (in layers of vertex neighbours) of each local cell to `seeds`,
# up to `max_layer`; -1 for cells further away
function cell_layers(lmodel,seeds,max_layer)
  topo = Gridap.Geometry.get_grid_topology(lmodel)
  Dc   = num_cell_dims(lmodel)
  cell_to_vertices = get_faces(topo,Dc,0)
  vertex_to_cells  = get_faces(topo,0,Dc)
  layer = fill(-1,num_cells(lmodel))
  layer[seeds] .= 0
  front = collect(Int,seeds)
  for k in 1:max_layer
    next = Int[]
    for cell in front, vertex in cell_to_vertices[cell], ncell in vertex_to_cells[vertex]
      if layer[ncell] == -1
        layer[ncell] = k
        push!(next,ncell)
      end
    end
    front = next
  end
  layer
end

# Non-conforming glue of `portion`, a DiscreteModelPortion of a local model
# with glue `ncglue`. Gridap numbers the d-faces (d < Dc) of the portion in
# increasing order of their parent ids, so regular faces still come first.
# Hanging faces stay hanging. If the owner cell of a hanging face is not in the
# portion, it becomes -1, as in models built by p4est for hanging faces whose
# owner cell is not local.
function restrict_non_conforming_glue(pXest_refinement_rule,ncglue,portion)
  Dc   = num_cell_dims(portion)
  topo = Gridap.Geometry.get_grid_topology(portion)
  cell_to_old_cell = Gridap.Geometry.get_cell_to_parent_cell(portion)
  old_cell_to_cell = zeros(Int,num_cells(Gridap.Geometry.get_parent_model(portion)))
  old_cell_to_cell[cell_to_old_cell] .= eachindex(cell_to_old_cell)

  num_regular_faces      = Vector{Int}(undef,Dc)
  num_hanging_faces      = Vector{Int}(undef,Dc)
  hanging_faces_glue     = Vector{Vector{Tuple{Int,Int,Int}}}(undef,Dc)
  hanging_faces_to_cell  = Vector{Vector{Int}}(undef,Dc)
  hanging_faces_to_lface = Vector{Vector{Int}}(undef,Dc)
  owner_faces_pindex     = Vector{Vector{Int}}(undef,Dc-1)
  owner_faces_lids       = Vector{Dict{Int,Tuple{Int,Int,Int}}}(undef,Dc-1)
  cell_faces             = Vector{Gridap.Arrays.Table}(undef,Dc)
  for d in 0:Dc-1
    face_to_old_face = Gridap.Geometry.get_face_to_parent_face(portion,d)
    num_regular_faces_old = ncglue.num_regular_faces[d+1]
    @assert issorted(face_to_old_face)
    num_regular_faces[d+1] = count(<=(num_regular_faces_old),face_to_old_face)
    num_hanging_faces[d+1] = length(face_to_old_face) - num_regular_faces[d+1]
    hanging_faces_glue[d+1] = map(face_to_old_face[num_regular_faces[d+1]+1:end]) do old_face
      ocell, ocell_lface, subface = ncglue.hanging_faces_glue[d+1][old_face-num_regular_faces_old]
      new_ocell = (ocell > 0 && old_cell_to_cell[ocell] > 0) ? old_cell_to_cell[ocell] : -1
      (new_ocell, ocell_lface, subface)
    end
    cell_faces[d+1] = Gridap.Arrays.Table(get_faces(topo,Dc,d))
    hanging_faces_to_cell[d+1], hanging_faces_to_lface[d+1] =
      GridapP4est._generate_hanging_faces_to_cell_and_lface(num_regular_faces[d+1],
                                                            num_hanging_faces[d+1],
                                                            cell_faces[d+1])
  end
  GridapP4est._fill_owner_faces_pindex_and_lids!(Dc,pXest_refinement_rule,
                                                 owner_faces_pindex,owner_faces_lids,
                                                 num_hanging_faces,hanging_faces_glue,
                                                 hanging_faces_to_cell,hanging_faces_to_lface,
                                                 cell_faces)
  GridapP4est.NonConformingGlue(num_regular_faces,num_hanging_faces,hanging_faces_glue,
                                hanging_faces_to_cell,hanging_faces_to_lface,
                                owner_faces_pindex,owner_faces_lids)
end

function in_fe_space(;order=2)
  u(x) = x[1]^order + x[2]^order
  f(x) = -Δ(u)(x)
  ud(x) = u(x)
  return (u, f, ud)
end

function out_fe_space(;kx=1.0,ky=1.0)
  u(x) = sin(2*pi*kx*x[1]) * cos(2*pi*ky*x[2])
  f(x) = -Δ(u)(x)
  ud(x) = u(x)
  return (u, f, ud)
end

function solve_on_model(model,cutgeo,geo,geometry::UnfittedGeometry;
                        order=2::Integer,
                        problem::Tuple=in_fe_space(),
                        write_solution::Bool=false,
                        vtk_suffix::String="")
  vtk_name(name::String) = isempty(vtk_suffix) ? name : "$(name)_$(vtk_suffix)"

  # Generate aggregates and the aggregate-conforming model (with one ghost layer)
  cell_to_root, agg_cell_indices = generate_aggregates(model,cutgeo,geo)
  model, new_to_old_cell = aggregate_conforming_model(model,agg_cell_indices)

  # The local cells of the new model are a renumbered subset of the old ones:
  # move the aggregates (root gids, unchanged) to the new local numbering, and
  # cut the new model (the old cutgeo is defined on the old local models)
  cell_to_root = map(cell_to_root,new_to_old_cell) do cell_to_root,new_to_old_cell
    cell_to_root[new_to_old_cell]
  end
  cutgeo = cut(model,geo)

  # The cell ownership now comes from the model, so the standard constructors
  # give the aggregate-conforming triangulations
  Ω = Triangulation(model)
  Ωa = Triangulation(GridapDistributed.WithGhost(),cutgeo,ACTIVE_IN)
  Ωp = Triangulation(GridapDistributed.NoGhost(),cutgeo,PHYSICAL_IN)
  Γ = EmbeddedBoundary(GridapDistributed.NoGhost(),cutgeo)
  n_Γ = get_normal_vector(Γ)

  ocell_to_root = map(cell_to_root,map(own_to_local,partition(get_cell_gids(model)))) do agg,o_to_l
    map(Reindex(agg),o_to_l)
  end

  if write_solution
    writevtk(model,datadir(vtk_name("model")));
    writevtk(Ω,datadir(vtk_name("triangulation")));
    writevtk(Γ,datadir(vtk_name("boundary")));
    writevtk(Ωa,datadir(vtk_name("active")));
    writevtk(Ωp,datadir(vtk_name("physical")));
    writevtk(Ω,
             datadir(vtk_name("aggregates")),
             celldata=["aggregate"=>ocell_to_root]);
  end

  reffe = ReferenceFE(lagrangian,Float64,order)
  spaces = map(local_views(Ωa)) do trian
    if geometry isa Union{SquareWithCoarseCuts,SinusoidalBand}
      FESpace(trian,reffe;conformity=:H1,dirichlet_tags=[1,2,3,4,7,8])
    else
      FESpace(trian,reffe;conformity=:H1)
    end
  end

  # 1. Generate AMR constraints (on the aggregate-conforming triangulation)

  amr_sDOF_to_dof, amr_sDOF_to_dofs, amr_sDOF_to_coeffs = 
    generate_amr_constraints(Ωa,spaces,reffe);

  # 2. Generate AGG constraints

  # 2.1. Generate AGG constraints (properly)
  agg_sDOF_gids, new_n_fdofs, new_n_ddofs, agg_mfdof_gids, agg_mddof_gids, agg_mDOF_to_dof, agg_sDOF_to_dof, agg_sDOF_to_mdofs, agg_sDOF_to_coeffs = 
    generate_merged_constraints(Ωa, 
                                spaces, 
                                cell_to_root, 
                                amr_sDOF_to_dof, 
                                amr_sDOF_to_dofs, # TODO: dof numbered and GD might assume DOF numbered
                                amr_sDOF_to_coeffs);

  u, f, ud = problem

  V = MyFESpace(
    Ωa,reffe,spaces,new_n_fdofs,new_n_ddofs,agg_mfdof_gids,
    agg_mDOF_to_dof,agg_sDOF_to_dof,agg_sDOF_to_mdofs,agg_sDOF_to_coeffs # TODO: Dirichlet BCs
  )

  uᵢ = interpolate_everywhere(u,V)
  U = TrialFESpace(V,ud)

  D = 2
  degree = 2*order*D
  dΩ = Measure(Ωp,degree)
  dΓ = Measure(Γ,degree)

  cell_meas = map(get_cell_measure∘Triangulation,local_views(Ω))
  meas = map(minimum,cell_meas) |> PartitionedArrays.getany
  h = meas^(1/D)
  γd = 10.0

  a(u,v) =
    ∫( ∇(v)⋅∇(u) )dΩ +
    ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )dΓ

  l(v) =
    ∫( v*f )dΩ + ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud )dΓ

  op = AffineFEOperator(a,l,U,V)
  uₕ = solve(op)

  if write_solution
    writevtk(Ωa,datadir(vtk_name("check_active")),
             cellfields=["ui"=>uᵢ,"ei"=>u-uᵢ,"uh"=>uₕ,"eh"=>u-uₕ],nsubcells=2);
    writevtk(Ωp,datadir(vtk_name("check_physical")),
             cellfields=["ui"=>uᵢ,"ei"=>u-uᵢ,"uh"=>uₕ,"eh"=>u-uₕ],nsubcells=2);
  end

  e = u - uₕ
  el2 = sqrt(sum(∫( e*e )dΩ))
  eh1 = sqrt(sum(∫( e*e + ∇(e)⋅∇(e) )dΩ))

  @show el2, eh1, h

  return el2, eh1, h
end

function run_single_test(distribute,
                         geometry::UnfittedGeometry;
                         num_parts::Integer=1,
                         initial_uniform_refs::Integer=4,
                         order::Integer=2,
                         problem::Tuple=in_fe_space(),
                         write_solution::Bool=false,
                         num_ghost_layers::Integer=1)
  ranks = distribute(LinearIndices((num_parts,)))
  model, cutgeo, geo =
    generate_unfitted_model(ranks,
                            num_parts,
                            geometry,
                            initial_uniform_refs=initial_uniform_refs,
                            num_ghost_layers=num_ghost_layers)
  el2, eh1, h = solve_on_model(
    model,cutgeo,geo,geometry;
    order,problem,write_solution
  )
  return el2, eh1, h
end

function run_convergence_test(geometry::UnfittedGeometry;
                              initial_uniform_refs::Integer,
                              order::Integer,
                              num_refinements::Integer,
                              problem::Tuple=in_fe_space(),
                              write_solution::Bool=false,
                              num_ghost_layers::Integer=1)
  
  @assert num_refinements >= 1 "num_refinements must be >= 1"

  el2s = Float64[]
  eh1s = Float64[]
  hs = Float64[]

  for level in 1:num_refinements
    suffix = "level_$(level)"
    model, cutgeo, geo = with_mpi() do distribute 
      ranks = distribute(LinearIndices((1,)))
      generate_unfitted_model(ranks,1,geometry,
        initial_uniform_refs=initial_uniform_refs,
        post_uniform_refinements=level-1,
        num_ghost_layers=num_ghost_layers)
    end
    el2, eh1, h = solve_on_model(
      model,cutgeo,geo,geometry;
      order,problem,write_solution,
      vtk_suffix=suffix
    )
    push!(el2s,el2)
    push!(eh1s,eh1)
    push!(hs,h)
  end

  function slope(hs,errors)
    x = log10.(hs)
    y = log10.(errors)
    linreg = hcat(fill!(similar(x),1),x) \ y
    linreg[2]
  end

  l2_slope = slope(hs,el2s)
  h1_slope = slope(hs,eh1s)

  p = plot(
    hs,
    [el2s eh1s],
    xaxis=:log,
    yaxis=:log,
    label=["L2 (slope = $(round(l2_slope,digits=2)))" "H1 (slope = $(round(h1_slope,digits=2)))"],
    shape=:auto,
    xlabel="h",
    ylabel="error norm",
    title="Convergence test",
    legend=:bottomright,
  )

  geom_name = lowercase(string(nameof(typeof(geometry))))
  plot_path = datadir("convergence_$(geom_name)_k$(order).png")
  mkpath(dirname(plot_path))
  savefig(p, plot_path)

  return el2s, eh1s, hs
end

end # module