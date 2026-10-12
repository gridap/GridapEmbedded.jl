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
                                 post_uniform_refinements=0)
  coarse_model = CartesianDiscreteModel((-1,1,-1,1),(1,1))
  num_ghost_layers = 1
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

function aggregate_conforming_model(model::OctreeDistributedDiscreteModel,agg_cell_indices)
  aggdmodel = GridapDistributed.GenericDistributedDiscreteModel(
    local_views(model),PRange(agg_cell_indices)
  )
  return OctreeDistributedDiscreteModel(
    num_cell_dims(model), num_point_dims(model),
    model.parts,aggdmodel,model.non_conforming_glue,
    model.coarse_model,
    model.ptr_pXest_connectivity,
    model.ptr_pXest,
    model.pXest_type,
    model.pXest_refinement_rule_type,
    model.connectivity_ref; num_ghost_layers = model.num_ghost_layers
  )
end

# Triangulation whose cell ownership is given by the aggregate-conforming cell indices
function aggregate_conforming_trian(TT,portion,model,cutgeo,agg_cell_indices,args...)
  trians = map(local_views(cutgeo),agg_cell_indices) do cutgeo,agg_cell_indices
    TT(portion,agg_cell_indices,cutgeo,args...)
  end
  return GridapDistributed.DistributedTriangulation(trians,model)
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

  if write_solution
    writevtk(Triangulation(model),datadir(vtk_name("old_model")));
  end

  # Generate aggregates and the aggregate-conforming triangulations
  cell_to_root, agg_cell_indices = generate_aggregates(model,cutgeo,geo)
  model = aggregate_conforming_model(model,agg_cell_indices)
  Ω = Triangulation(model)
  Ωa = aggregate_conforming_trian(Triangulation,GridapDistributed.WithGhost(),model,cutgeo,agg_cell_indices,ACTIVE_IN)
  Ωp = aggregate_conforming_trian(Triangulation,GridapDistributed.NoGhost(),model,cutgeo,agg_cell_indices,PHYSICAL_IN)
  Γ = aggregate_conforming_trian(EmbeddedBoundary,GridapDistributed.NoGhost(),model,cutgeo,agg_cell_indices)
  n_Γ = get_normal_vector(Γ)

  ocell_to_root = map(cell_to_root,map(own_to_local,agg_cell_indices)) do agg,o_to_l
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
                         write_solution::Bool=false)
  ranks = distribute(LinearIndices((num_parts,)))
  model, cutgeo, geo =
    generate_unfitted_model(ranks,
                            num_parts,
                            geometry,
                            initial_uniform_refs=initial_uniform_refs)
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
                              write_solution::Bool=false)
  
  @assert num_refinements >= 1 "num_refinements must be >= 1"

  el2s = Float64[]
  eh1s = Float64[]
  hs = Float64[]

  for level in 1:num_refinements
    suffix = "level_$(level)"
    model, cutgeo, geo = with_mpi() do distribute 
      generate_unfitted_model(distribute,geometry,
        initial_uniform_refs=initial_uniform_refs,
        post_uniform_refinements=level-1)
    end
    _, _, _, _, el2, eh1, h = solve_on_model(
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