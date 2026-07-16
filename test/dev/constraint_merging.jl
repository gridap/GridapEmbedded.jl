module DistributedAggregationP4estMeshes

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

function square()
  quadrilateral(;x0=Point(-0.9,-0.9),
                 d1=VectorValue(1.8,0.0),
                 d2=VectorValue(0.0,1.8))
end

function square_with_circular_hole()
  geo1 = quadrilateral(;x0=Point(-0.9,-0.9),
                        d1=VectorValue(1.8,0.0),
                        d2=VectorValue(0.0,1.8))
  geo2 = ! disk(0.4)
  geo = intersect(geo1,geo2)
  return geo
end

function flower(;x₀=Point(0.0,0.0),R₀=0.6,m=0.6,ω=10.0)
  name="flower"
  function flowerfun(x)
    _flower(x,x₀,R₀,m,ω)
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

function generate_unfitted_model(distribute,np,generate_geometry)
  ranks = distribute(LinearIndices((np,)))
  coarse_model = CartesianDiscreteModel((-1,1,-1,1),(1,1))
  num_uniform_refinements = 6
  num_ghost_layers = 1

  dmodel = OctreeDistributedDiscreteModel(
    ranks,coarse_model,num_uniform_refinements;num_ghost_layers=num_ghost_layers
  )

  geo = generate_geometry()
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

  for i in 1:num_uniform_refinements-1
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
  if np > 1
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

model, cutgeo, geo = with_mpi() do distribute
  generate_unfitted_model(distribute,1,flower)
end

Ω = Triangulation(model)
Γ = EmbeddedBoundary(cutgeo)
n_Γ = get_normal_vector(Γ)
Ωa = Triangulation(cutgeo,ACTIVE_IN)
Ωp = Triangulation(cutgeo,PHYSICAL_IN)

writevtk(Γ,datadir("boundary"));
writevtk(Ω,datadir("model"));
writevtk(Ωa,datadir("active"));
writevtk(Ωp,datadir("physical"));

cell_gids = partition(get_cell_gids(model))
cell_to_root = generate_aggregates(model,cutgeo,geo)

writevtk(Triangulation(model),datadir("aggregates"),celldata=["aggregate"=>cell_to_root]);

order = 2
reffe = ReferenceFE(lagrangian,Float64,order)
spaces = map(local_views(Ωa)) do trian
  FESpace(trian,reffe;conformity=:H1)
end

amr_sDOF_to_dof, amr_sDOF_to_dofs, amr_sDOF_to_coeffs = generate_amr_constraints(model,Ωa,spaces,reffe);

agg_sDOF_gids, agg_mfdof_gids, agg_mddof_gids, agg_mDOF_to_dof, agg_sDOF_to_dof, agg_sDOF_to_mdofs, agg_sDOF_to_coeffs = 
  generate_agfem_constraints(Ωa,spaces,cell_to_root); 

# Only works because we dont have bcs or exterior MDOFs, otherwise a bit more work is needed
agg_sDOF_to_dofs = map(agg_mDOF_to_dof, agg_sDOF_to_mdofs) do mDOF_to_dof, sDOF_to_mdofs
  data = mDOF_to_dof[sDOF_to_mdofs.data]
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
  is_master_of_well_posed_hanging = falses(ndofs)
  for sDOF in amr_only_sDOF
    masters = amr_sDOF_to_dofs[sDOF]
    is_master_of_well_posed_hanging[masters] .= true
  end
  keep_agg_constraint = falses(length(agg_sDOF_to_dof))
  for (sDOF,dof) in enumerate(agg_sDOF_to_dof)
    keep_agg_constraint[sDOF] = !is_master_of_well_posed_hanging[dof]
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

u(x) = x[1]^2 + x[2]^2
f(x) = -Δ(u)(x)
ud(x) = u(x)

V = MyFESpace(Ωa,reffe,spaces,
              sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs)
uᵢ = interpolate_everywhere(u,V)
U = TrialFESpace(V)

D = 2
degree = 2*order*D
dΩ = Measure(Ωp,degree)
dΓ = Measure(Γ,degree)

cell_meas = map(get_cell_measure∘Triangulation,local_views(Ω))
meas = map(minimum,cell_meas) |> PartitionedArrays.getany
const h = meas^(1/D)
const γd = 10.0

a(u,v) =
  ∫( ∇(v)⋅∇(u) )dΩ +
  ∫( (γd/h)*v*u  - v*(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))*u )dΓ

l(v) =
  ∫( v*f )dΩ + ∫( (γd/h)*v*ud - (n_Γ⋅∇(v))*ud )dΓ

op = AffineFEOperator(a,l,U,V)
uₕ = solve(op)

writevtk(Ωa,datadir("check_active"),
         cellfields=["ui"=>uᵢ,"ei"=>u-uᵢ,"uh"=>uₕ,"eh"=>u-uₕ]);
writevtk(Ωp,datadir("check_physical"),
         cellfields=["ui"=>uᵢ,"ei"=>u-uᵢ,"uh"=>uₕ,"eh"=>u-uₕ]);

end # module