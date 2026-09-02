module DistributedAgFEM4P4estMeshes

using Gridap
using Gridap.Helpers
using GridapEmbedded
using GridapEmbedded.Interfaces: CUT
using GridapEmbedded.LevelSetCutters: Leaf
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

# This respects the input triangulation, without adding back all the ghosts
function MyFESpace(
  trian::GridapDistributed.DistributedTriangulation,reffe,spaces,
  sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs;split_own_and_ghost=false,constraint=nothing,kwargs...
)
  myfespaces = map(sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs,spaces) do sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs,spaces
    FESpaceWithLinearConstraints(sDOF_to_dof,
                                 sDOF_to_dofs,
                                 sDOF_to_coeffs,
                                 spaces)
  end
  gids = GridapDistributed.generate_gids(trian,myfespaces)
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
  quadrilateral(;x0=Point(-1.1,-0.88),
                 d1=VectorValue(2.2,0.0),
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

function generate_unfitted_model(distribute,
                                 geometry::UnfittedGeometry;
                                 num_parts=1,
                                 initial_uniform_refs=4,
                                 post_uniform_refinements=0)
  ranks = distribute(LinearIndices((num_parts,)))
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
  if num_parts > 1
    fmodel, _ = GridapDistributed.redistribute(fmodel)
  end
  fmodel = Gridap.Adaptivity.get_model(fmodel)
  cutgeo = cut(fmodel,geo)
  return fmodel, cutgeo, geo
end

function generate_amr_constraints(model, trian, spaces, reffe)
  cell_gids = get_cell_gids(model)
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
  ncgt = NonConformingGridTopology(fmodel)
  strategy = AggregateCutCellsByThreshold(1.0)
  _, lcell_to_root, _ = map(
    local_views(cutgeo),partition(cell_gids),ncgt
  ) do cutgeo,cell_indices,ncgt
    lid_to_gid = local_to_global(cell_indices)
    aggregate(strategy,cutgeo,geo,lid_to_gid,IN,grid_topology=ncgt)
  end |> tuple_of_arrays
  return lcell_to_root
end

function generate_agfem_constraints(trian, spaces, bgcell_to_bgroot)
  shfns = DistributedCellField(map(get_fe_basis, spaces), trian)
  dofs = DistributedCellDof(map(get_fe_dof_basis, spaces), trian)
  bgcell_to_gcell = map(local_views(trian)) do trian
    glue = get_glue(trian,Val(num_cell_dims(trian)))
    glue.mface_to_tface
  end
  sDOF_gids, mfdof_gids, mddof_gids, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs = 
  GridapEmbedded.Distributed.generate_aggregated_space_constraints(
    trian, spaces, bgcell_to_bgroot, shfns, dofs, bgcell_to_gcell
  )
  return sDOF_gids, mfdof_gids, mddof_gids, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs
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
  Ω = Triangulation(model)
  Γ = EmbeddedBoundary(cutgeo)
  n_Γ = get_normal_vector(Γ)
  Ωa = Triangulation(cutgeo,ACTIVE_IN)
  Ωp = Triangulation(cutgeo,PHYSICAL_IN)

  if write_solution
    writevtk(Γ,datadir(vtk_name("boundary")));
    writevtk(Ω,datadir(vtk_name("model")));
    writevtk(Ωa,datadir(vtk_name("active")));
    writevtk(Ωp,datadir(vtk_name("physical")));
  end

  cell_to_root = generate_aggregates(model,cutgeo,geo)

  if write_solution
    writevtk(Triangulation(model),datadir(vtk_name("aggregates")),celldata=["aggregate"=>cell_to_root]);
  end

  reffe = ReferenceFE(lagrangian,Float64,order)
  spaces = map(local_views(Ωa)) do trian
    if geometry isa Union{SquareWithCoarseCuts,SinusoidalBand}
      FESpace(trian,reffe;conformity=:H1,dirichlet_tags=[1,2,3,4,7,8])
    else
      FESpace(trian,reffe;conformity=:H1)
    end
  end

  amr_sDOF_to_dof, amr_sDOF_to_dofs, amr_sDOF_to_coeffs = generate_amr_constraints(model,Ωa,spaces,reffe);

  agg_sDOF_gids, agg_mfdof_gids, agg_mddof_gids, agg_mDOF_to_dof, agg_sDOF_to_dof, agg_sDOF_to_mdofs, agg_sDOF_to_coeffs = 
    generate_agfem_constraints(Ωa,spaces,cell_to_root);

  # Deal with Dirichlet BCs
  agg_sDOF_to_dofs = map(agg_mDOF_to_dof, agg_sDOF_to_mdofs) do mDOF_to_dof, sDOF_to_mdofs
    T = eltype(sDOF_to_mdofs.data)
    data = Vector{T}(undef, length(sDOF_to_mdofs.data))
    @inbounds for i in eachindex(sDOF_to_mdofs.data)
      mdof = sDOF_to_mdofs.data[i]
      data[i] = ifelse(mdof > zero(T), mDOF_to_dof[mdof], mdof)
    end
    return Table(data, sDOF_to_mdofs.ptrs)
  end

  # Remove ill-posed constraints that are also hanging
  agg_sDOF_keep = map(
    spaces, amr_sDOF_to_dof, amr_sDOF_to_dofs, agg_sDOF_to_dof, agg_sDOF_to_dofs,
  ) do space, amr_sDOF_to_dof, amr_sDOF_to_dofs, agg_sDOF_to_dof, agg_sDOF_to_dofs
    ndofs = num_free_dofs(space)
    is_agg_dof = falses(ndofs)
    is_agg_dof[agg_sDOF_to_dof] .= true
    amr_only_sDOF = findall(dof->!is_agg_dof[dof],amr_sDOF_to_dof)
    is_mfdof_of_well_posed_hanging = falses(ndofs)
    for sDOF in amr_only_sDOF
      masters = amr_sDOF_to_dofs[sDOF]
      fmasters = masters[masters .> 0]
      is_mfdof_of_well_posed_hanging[fmasters] .= true
    end
    keep_agg_constraint = falses(length(agg_sDOF_to_dof))
    for (sDOF,dof) in enumerate(agg_sDOF_to_dof)
      keep_agg_constraint[sDOF] = !is_mfdof_of_well_posed_hanging[dof]
    end
    display("Removing $(count(!,keep_agg_constraint)) agg constraints:")
    display(findall(!,keep_agg_constraint))
    return findall(keep_agg_constraint)
  end;
  agg_sDOF_to_dof, agg_sDOF_to_dofs, agg_sDOF_to_coeffs = map(
    agg_sDOF_keep, agg_sDOF_to_dof, agg_sDOF_to_dofs, agg_sDOF_to_coeffs
  ) do agg_sDOF_keep, agg_sDOF_to_dof, agg_sDOF_to_dofs, agg_sDOF_to_coeffs
    return agg_sDOF_to_dof[agg_sDOF_keep], agg_sDOF_to_dofs[agg_sDOF_keep], agg_sDOF_to_coeffs[agg_sDOF_keep]
  end |> tuple_of_arrays;

  # Merge and close the constraint tables, preferring the AMR constraints in case of conflict 
  # for ill-posed dofs that are also hanging
  sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = map(
    spaces, agg_sDOF_to_dof, agg_sDOF_to_dofs, agg_sDOF_to_coeffs, amr_sDOF_to_dof, amr_sDOF_to_dofs, amr_sDOF_to_coeffs
  ) do space, agg_sDOF_to_dof, agg_sDOF_to_dofs, agg_sDOF_to_coeffs, amr_sDOF_to_dof, amr_sDOF_to_dofs, amr_sDOF_to_coeffs
    on_conflict(dof, dofs1, coeffs1, offset1, dofs2, coeffs2, offset2) = (dofs2, coeffs2, offset2)
    sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, sDOF_to_offsets = FESpaces.merge_slave_constraint_tables(
      space, 
      agg_sDOF_to_dof, agg_sDOF_to_dofs, agg_sDOF_to_coeffs, 
      amr_sDOF_to_dof, amr_sDOF_to_dofs, amr_sDOF_to_coeffs;
      on_conflict
    )
    sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, sDOF_to_offsets = FESpaces.close_slave_constraint_tables(
      space, sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, sDOF_to_offsets
    )
  end |> tuple_of_arrays;

  u, f, ud = problem

  V = MyFESpace(Ωa,reffe,spaces,
                sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs)
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
             cellfields=["ui"=>uᵢ,"ei"=>u-uᵢ,"uh"=>uₕ,"eh"=>u-uₕ]);
    writevtk(Ωp,datadir(vtk_name("check_physical")),
             cellfields=["ui"=>uᵢ,"ei"=>u-uᵢ,"uh"=>uₕ,"eh"=>u-uₕ]);
  end

  e = u - uₕ
  el2 = sqrt(sum(∫( e*e )dΩ))
  eh1 = sqrt(sum(∫( e*e + ∇(e)⋅∇(e) )dΩ))

  return model, cutgeo, geo, uₕ, el2, eh1, h
end

function run_single_test(geometry::UnfittedGeometry;
                         initial_uniform_refs::Integer,
                         order::Integer,
                         problem::Tuple=in_fe_space(),
                         write_solution::Bool=false)
  model, cutgeo, geo = with_mpi() do distribute
    generate_unfitted_model(distribute,geometry,initial_uniform_refs=initial_uniform_refs)
  end
  model, cutgeo, geo, uₕ, _, _, _ = solve_on_model(
    model,cutgeo,geo,geometry;
    order,
    problem,
    write_solution
  )
  return model, cutgeo, geo, uₕ
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
      order,
      problem,
      write_solution,
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

# Square:                 begin with initial_uniform_refs = 3
# SquareWithCircularHole: begin with initial_uniform_refs = 3
# Flower:                 begin with initial_uniform_refs = 6
# SinusoidalBand:         begin with initial_uniform_refs = 5
el2s, eh1s, hs = run_convergence_test(SquareWithCircularHole(),
                     initial_uniform_refs=3,
                     order=2,
                     num_refinements=4,
                     problem=out_fe_space())

end # module