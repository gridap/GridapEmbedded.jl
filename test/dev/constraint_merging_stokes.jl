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

# REFERENCES
#
# 1. p4est: SCALABLE ALGORITHMS FOR PARALLEL ADAPTIVE MESH REFINEMENT ON FORESTS OF OCTREES
# 2. A GENERIC FINITE ELEMENT FRAMEWORK ON PARALLEL TREE-BASED ADAPTIVE MESHES https://arxiv.org/pdf/1907.03709
# 3. THE AGGREGATED UNFITTED FINITE ELEMENT METHOD ON PARALLEL TREE-BASED ADAPTIVE MESHES https://arxiv.org/pdf/2006.05373

# BRANCHES
# 1. GridapP4est.jl: https://github.com/gridap/GridapP4est.jl/tree/fe_space_on_triangulation
# 2. Gridap.jl: https://github.com/gridap/Gridap.jl/tree/constraints
# 3. GridapDistributed.jl: https://github.com/gridap/GridapDistributed.jl/tree/agfem
# 4. GridapEmbedded.jl: https://github.com/gridap/GridapEmbedded.jl/tree/distributed_aggregate_p4est_meshes

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
    model.owns_ptr_pXest_connectivity,
    model.gc_ref; num_ghost_layers = model.num_ghost_layers
  )
end

include("NonConformingGridTopologies.jl")
using .NonConformingGridTopologies

function generate_unfitted_model(distribute,geo;
                                 num_parts=1,initial_uniform_refs=4)
  ranks = distribute(LinearIndices((num_parts,)))
  coarse_model = CartesianDiscreteModel((-1,1,0,1),(2,1))
  num_ghost_layers = 1
  dmodel = OctreeDistributedDiscreteModel(
    ranks,coarse_model,initial_uniform_refs;num_ghost_layers=num_ghost_layers
  )
  # From here on we illustrate how to refine and coarsen the model based on the cut geometry
  # Beware that cell_to_inoutcut and fmodel_refine_coarsen_flags are distributed MPI arrays
  # https://partitionedarrays.github.io/PartitionedArrays.jl/stable/usage/#Basic-usage
  # See also tests from PartitionedArrays.jl and GridapDistributed.jl for more examples
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
  # Exercice: Plot the model at this stage to see the effect of the refinement and coarsening flags
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
    aggregate(strategy,cutgeo,geo,lid_to_gid,OUT,grid_topology=ncgt)
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

# My hope is that you can use this function as a blackbox
function generate_constrained_fe_space(model, Ωa, spaces, reffe, cell_to_root)
  amr_sDOF_to_dof, amr_sDOF_to_dofs, amr_sDOF_to_coeffs = generate_amr_constraints(model,Ωa,spaces,reffe);
  _, _, _, agg_mDOF_to_dof, agg_sDOF_to_dof, agg_sDOF_to_mdofs, agg_sDOF_to_coeffs =
    generate_agfem_constraints(Ωa,spaces,cell_to_root);

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

  V = MyFESpace(Ωa,reffe,spaces,
                sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs)
  return V
end

function in_fe_space_if_u_order_two()
  u(x) = VectorValue( x[1]*x[1], x[2] )
  p(x) = x[1]
  f(x) = -Δ(u)(x) + ∇(p)(x)
  g(x) = (∇⋅u)(x)
  ∇u(x) = ∇(u)(x)
  return (u, p, f, g, ∇u)
end

function out_fe_space()
  @notimplemented "Out FE space example is not implemented yet"
end

function solve_on_model(model,cutgeo,geo;
                        order=2::Integer,
                        problem::Tuple=in_fe_space_if_u_order_two(),
                        write_solution::Bool=false)

  # Update à la algoim way                     
  Ω = Triangulation(model)
  Γ = EmbeddedBoundary(cutgeo)
  nΓ = get_normal_vector(Γ)
  Ωa = Triangulation(cutgeo,ACTIVE_OUT)
  Ωp = Triangulation(cutgeo,PHYSICAL_OUT)

  if write_solution
    writevtk(Γ,datadir("boundary"));
    writevtk(Ω,datadir("model"));
    writevtk(Ωa,datadir("active"));
    writevtk(Ωp,datadir("physical"));
  end

  cell_to_root = generate_aggregates(model,cutgeo,geo)

  if write_solution
    writevtk(Triangulation(model),datadir("aggregates"),celldata=["aggregate"=>cell_to_root]);
  end

  D = 2
  
  reffeᵘ = ReferenceFE(lagrangian,VectorValue{D,Float64},order)
  spacesᵘ = map(local_views(Ωa)) do trian
    TestFESpace(trian,reffeᵘ;conformity=:H1,dirichlet_tags="boundary")
      # dirichlet_tags=[1,2,5, 7,8, 3,4,6],
      # dirichlet_masks=[(false,true),(false,true),(false,true),
      #                  (false,true),(false,true),
      #                  (false,true),(false,true),(false,true)])
  end
  V = generate_constrained_fe_space(model, Ωa, spacesᵘ, reffeᵘ, cell_to_root)

  reffeᵖ = ReferenceFE(lagrangian,Float64,order-1,space=:P)
  spacesᵖ = map(local_views(Ωa)) do trian
    TestFESpace(trian,reffeᵖ)
  end
  Q = generate_constrained_fe_space(model, Ωa, spacesᵖ, reffeᵖ, cell_to_root)

  K = ConstantFESpace(model)

  u, p, f, g, ∇u = problem

  U = TrialFESpace(V,u)
  P = TrialFESpace(Q)
  L = TrialFESpace(K)

  Y = MultiFieldFESpace([V,Q,K])
  X = MultiFieldFESpace([U,P,L])

  # First sanity check: interpolate a FE function that can be exactly represented by the FE space that combines hanging and aggregate constraints
  uᵢ, pᵢ, _ = interpolate_everywhere([u,p,0.0],X)

  D = 2
  # Update with algoim
  degree = 2*order*D
  dΩ = Measure(Ωp,degree)
  dΓ = Measure(Γ,degree)

  cell_meas = map(get_cell_measure∘Triangulation,local_views(Ω))
  meas = map(minimum,cell_meas) |> PartitionedArrays.getany
  h = meas^(1/D)
  γ = 10.0

  # [WARNING] Sign of normal is flipped because OUT domain

  α(u,v) = ∫( ∇(u)⊙∇(v) )dΩ
  β(v,q) = ∫( q*(∇⋅v) )dΩ

  γˡ(u,v,p,q) = ∫( (γ/h)*(v⋅u) +
    v⋅(nΓ⋅∇(u)) + (nΓ⋅∇(v))⋅u - (p*nΓ)⋅v - (q*nΓ)⋅u )dΓ
  γʳ(v,q) = ∫( (γ/h)*(v⋅u) + (nΓ⋅∇(v))⋅u - (q*nΓ)⋅u )dΓ

  r(p,ℓ) = ∫( p*ℓ )dΩ

  a((u,p,l),(v,q,ℓ)) = 
    α(u,v) - β(v,p) - β(u,q) + 
    γˡ(u,v,p,q) + r(p,ℓ) + r(q,l)
  b((v,q,l)) = γʳ(v,q) + ∫( v⋅f - q*g )dΩ

  op = AffineFEOperator(a,b,X,Y)
  uₕ, pₕ, _ = solve(op)

  # Second sanity check: If the first sanity check passes, but the solution to the problem is incorrect, the most typical cause is a problem in the definition of the weak form, e.g. a missing term or a wrong sign.

  if write_solution
    writevtk(Ωa,datadir("check_active"),
            cellfields=["ui"=>uᵢ,"eui"=>u-uᵢ,"p"=>pᵢ,"epi"=>p-pᵢ,
                        "uh"=>uₕ,"euh"=>u-uₕ,"ph"=>pₕ,"eph"=>p-pₕ]);
    writevtk(Ωp,datadir("check_physical"),
            cellfields=["ui"=>uᵢ,"eui"=>u-uᵢ,"p"=>pᵢ,"epi"=>p-pᵢ,
                        "uh"=>uₕ,"euh"=>u-uₕ,"ph"=>pₕ,"eph"=>p-pₕ]);
  end

end

# MAIN PROGRAMME STARTS HERE

geo = disk(0.45)
initial_uniform_refs = 4
order = 2
problem = in_fe_space_if_u_order_two()
write_solution = true

# Update generate_unfitted_model to
# generate model with an algoim LS
# HINT: Use s_cell_quad,is_c₋ = CellQuadratureAndActiveMask(bgmodel,squad)
# The model with cutgeo is generate by coarsening a maximumly refined model, so that the cut geometry is well resolved
# If you start with an analytical LS of disk, you can start with a coarse model (initial_uniform_refs=1) and refine on cut cells until achieving the desired mesh size on the cuts
model, cutgeo, _ = with_mpi() do distribute
  generate_unfitted_model(distribute,geo,initial_uniform_refs=initial_uniform_refs)
end

# Although it does not run on a distributed environment,
# the input arrays are distributed, so all the derived
# types from them are distributed 
solve_on_model(
  model, cutgeo, geo;
  order,
  problem,
  write_solution
)

end # module