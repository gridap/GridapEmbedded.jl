module DistributedAggregationP4estMeshes

  using Gridap
  using GridapEmbedded
  using GridapEmbedded.Interfaces: CUT
  using GridapDistributed
  using PartitionedArrays
  using MPI

  using Gridap.Geometry: get_faces

  using P4est_wrapper
  using GridapP4est
  using GridapP4est: generate_local_fe_spaces_and_constraints

  include("NonConformingGridTopologies.jl")
  using .NonConformingGridTopologies

  function run(distribute)

    ranks = distribute(LinearIndices((MPI.Comm_size(MPI.COMM_WORLD),)))
    coarse_model = CartesianDiscreteModel((-1,-0.5,-1,-0.5),(1,1))
    num_uniform_refinements = 1
    num_ghost_layers = 1

    dmodel = OctreeDistributedDiscreteModel(ranks,
                                            coarse_model,
                                            num_uniform_refinements;
                                            num_ghost_layers=num_ghost_layers)
    
    geo1 = quadrilateral(;x0=Point(-0.9,-0.9),
                         d1=VectorValue(1.8,0.0),
                         d2=VectorValue(0.0,1.8))
    geo2 = ! disk(0.4)
    geo = geo1 # intersect(geo1,geo2)

    cutgeo = cut(dmodel, geo)
    cell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geo)
    fmodel_refine_coarsen_flags = 
      map(ranks,
          partition(get_cell_gids(dmodel.dmodel)),
          cell_to_inoutcut) do rank,indices,cell_to_inoutcut
        flags = zeros(Int,length(indices))
        flags .= nothing_flag
        toref = findall(c->c==CUT,cell_to_inoutcut)
        flags[toref] .= refine_flag
        flags
    end
    fmodel,_ = Gridap.Adaptivity.adapt(dmodel,fmodel_refine_coarsen_flags);

    for i in 1:1 # 5
      cutgeo = cut(fmodel, geo)
      cell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geo)
      fmodel_refine_coarsen_flags = 
        map(partition(get_cell_gids(fmodel)),
            cell_to_inoutcut) do indices,cell_to_inoutcut
          flags = zeros(Int,length(indices))
          flags .= nothing_flag
          toref = findall(c->c==CUT,cell_to_inoutcut)
          flags[toref] .= refine_flag
          flags
      end
      fmodel,_ = Gridap.Adaptivity.adapt(fmodel,fmodel_refine_coarsen_flags);
    end
    fmodel, _ = GridapDistributed.redistribute(fmodel)

    cutgeo = cut(fmodel, geo)
    Γ = EmbeddedBoundary(cutgeo)
    Ωᵃ = Triangulation(cutgeo,ACTIVE_IN)
    Ωᵖ = Triangulation(cutgeo,PHYSICAL_IN)

    writevtk(Γ,"data/quad_bnd");
    # The next line fails unless I comment out the assertion 
    # in the constructor of AppendedTriangulations.jl:
    # @assert get_background_model(a) === get_background_model(b)
    writevtk(Ωᵃ,"data/quad_act");
    # writevtk(Ωᵖ,"data/quad_phys");

    cell_gids = get_cell_gids(fmodel)
    cell_indices = partition(cell_gids)

    ncgt = NonConformingGridTopology(fmodel)
    strategy = AggregateCutCellsByThreshold(1.0)
    lcell_to_lroot, lcell_to_root, lcell_to_value =
      map(local_views(cutgeo),cell_indices,ncgt) do cutgeo,cell_indices,ncgt
        lid_to_gid = local_to_global(cell_indices)
        aggregate(strategy,cutgeo,geo,lid_to_gid,IN,grid_topology=ncgt)
      end |> tuple_of_arrays

    # Output for verification of lcell_to_root map
    ocell_to_root = map(lcell_to_root,own_to_local(cell_gids)) do agg,o_to_l
      map(Reindex(agg),o_to_l)
    end

    writevtk(Γ,"data/quad_bnd");
    writevtk(Ωᵖ,"data/quad_phys");
    writevtk(Triangulation(fmodel),"data/quad_aggregates",
             celldata = ["aggregate"=>ocell_to_root]);

    order = 2
    u(x) = x[1]^2 + x[2]^2
    reffe = ReferenceFE(lagrangian,Float64,order)
    
    # 0. FE space without resolved hanging dof constraints (i.e. with ill-posed free dofs)
    _, sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, _, _, _, _ =
        generate_local_fe_spaces_and_constraints(
          Ωᵃ,reffe;conformity=:H1,dirichlet_tags="boundary")
    # 1. FE space with resolved hanging dof constraints
    Vₕ = FESpace(Ωᵃ,reffe,conformity=:H1;dirichlet_tags="boundary")
    Uₕ = TrialFESpace(Vₕ,u)
    # A simple sanity check
    uₕ = interpolate_everywhere(u,Uₕ)
    writevtk(Ωᵖ,"data/quad_interpolated",cellfields=["uh"=>uₕ]);

    map(local_views(Vₕ),
        lcell_to_root,
        sDOF_to_dof,
        sDOF_to_dofs,
        sDOF_to_coeffs) do Vₕ,aggregates,sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs
      hagfemspace = hAgFEMSpace(Vₕ,aggregates,sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs)
      uₕ = interpolate_everywhere(u,hagfemspace)
      writevtk(Ωᵖ,"data/quad_check",cellfields=["uh"=>uₕ,"eh"=>u-uₕ]);
    end

    return nothing
  end

  using Gridap
  using Gridap.Helpers
  using Gridap.Arrays
  using Gridap.Geometry
  using Gridap.Geometry: get_cell_to_parent_cell
  using Gridap.ReferenceFEs
  using Gridap.CellData
  using Gridap.FESpaces

  function hAgFEMSpace(
    f::SingleFieldFESpace,
    bgcell_to_bgcellin::AbstractVector,
    sDOF_to_dof::AbstractVector,
    sDOF_to_dofs::AbstractVector,
    sDOF_to_coeffs::AbstractVector,
    g::SingleFieldFESpace=f,
    args...)

    @assert get_triangulation(f) === get_triangulation(g)
    hAgFEMSpace(f,bgcell_to_bgcellin,get_fe_basis(g),get_fe_dof_basis(g),sDOF_to_dof,sDOF_to_dofs,sDOF_to_coeffs,args...)
  end

  # Note: cell is in fact bgcell in this function since f will usually be an ExtendedFESpace
  function hAgFEMSpace(
    f::SingleFieldFESpace,
    bgcell_to_bgcellin::AbstractVector,
    shfns_g::CellField,
    dofs_g::CellDof,
    sDOF_to_dof::AbstractVector,
    sDOF_to_dofs::AbstractVector,
    sDOF_to_coeffs::AbstractVector,
    bgcell_to_gcell::AbstractVector=1:length(bgcell_to_bgcellin))

    # Triangulation made of active cells
    trian_a = get_triangulation(f)

    # Build root cell map (i.e. aggregates) in terms of active cell ids
    D = num_cell_dims(trian_a)
    glue = get_glue(trian_a,Val(D))
    acell_to_bgcell = glue.tface_to_mface
    bgcell_to_acell = glue.mface_to_tface
    acell_to_bgcellin = collect(lazy_map(Reindex(bgcell_to_bgcellin),acell_to_bgcell))
    acell_to_acellin = collect(lazy_map(Reindex(bgcell_to_acell),acell_to_bgcellin))
    acell_to_gcell = lazy_map(Reindex(bgcell_to_gcell),acell_to_bgcell)

    # Build shape funs of g by replacing local funs in cut cells by the ones at the root
    # This needs to be done with shape functions in the physical domain
    # otherwise shape funs in cut and root cells are the same
    acell_phys_shapefuns_g = get_array(change_domain(shfns_g,PhysicalDomain()))
    acell_phys_root_shapefuns_g = lazy_map(Reindex(acell_phys_shapefuns_g),acell_to_acellin)
    root_shfns_g = GenericCellField(acell_phys_root_shapefuns_g,trian_a,PhysicalDomain())

    # Compute data needed to compute the constraints
    dofs_f = get_fe_dof_basis(f)
    shfns_f = get_fe_basis(f)
    acell_to_coeffs = dofs_f(root_shfns_g)
    acell_to_proj = dofs_g(shfns_f)
    acell_to_dof_ids = get_cell_dof_ids(f)
    # Need DoF numbering from original non-conforming space
    acell_to_nc_dof_ids = get_cell_dof_ids(f.space)
    # Need constraint map from root ldofs to lmdofs to resolve
    # the hanging dof constraints of ill-posed free dofs
    acellin_constraints = lazy_map(Reindex(get_cell_constraints(f)),acell_to_acellin)

    aggdof_to_dof, aggdof_to_dofs, aggdof_to_coeffs = _setup_hagfem_constraints(
      num_free_dofs(f.space),
      num_free_dofs(f),
      num_free_dofs(f.space)-num_free_dofs(f),
      acell_to_acellin,
      acell_to_dof_ids,
      acell_to_nc_dof_ids,
      acell_to_coeffs,
      acell_to_proj,
      acell_to_gcell,
      f.mDOF_to_DOF,
      f.DOF_to_mDOFs,
      f.DOF_to_coeffs,
      sDOF_to_dof,
      sDOF_to_dofs,
      sDOF_to_coeffs,
      acellin_constraints) # Extract constraint data

    # Throwing away initial FESpaceWithLinearConstraints 
    # where only the hanging dofs were constrained
    FESpaceWithLinearConstraints(aggdof_to_dof,aggdof_to_dofs,aggdof_to_coeffs,f.space)
  end

  using GridapEmbedded.AgFEM: _allocate_fdof_to_data
  using GridapEmbedded.AgFEM: _fill_aggdof_to_dofs_ptrs!
  using GridapEmbedded.AgFEM: _allocate_aggdof_to_data
  import GridapEmbedded.AgFEM: _fill_aggdof_to_dofs_data!
  import GridapEmbedded.AgFEM: _fill_aggdof_to_coeffs_data!

  function _setup_hagfem_constraints(
    n_fdofs,  # Number of free dofs in the non-conforming space
    n_fmdofs, # Number of free dofs in the conforming space
    n_hdofs,  # Number of hanging dofs
    acell_to_acellin,
    acell_to_dof_ids,
    acell_to_nc_dof_ids,
    acell_to_coeffs,
    acell_to_proj,
    acell_to_gcell,
    mDOF_to_DOF,
    DOF_to_mDOFs,
    DOF_to_coeffs,
    sDOF_to_dof,
    sDOF_to_dofs,
    sDOF_to_coeffs,
    acellin_constraints,
    acell_to_is_owned=fill(true,length(acell_to_acellin)))

    # REFERENCE: [R] https://arxiv.org/pdf/2006.05373

    # Refactor of FESpaceWithLinearConstraints:
    # https://github.com/gridap/Gridap.jl/blob/constraints/src/FESpaces/FESpacesWithLinearConstraints.jl

    # # #
    # 1. Identify ill-posed free dofs (fdof_to_is_agg), 
    #    their root cells (fdof_to_acell) and 
    #    the local dof number on the root cell (fdof_to_ldof)
    # # #

    ##### OLD

    # [!] fdofs follow the numbering of the conforming space
    fdof_to_is_agg, fdof_to_acell, fdof_to_ldof = 
      _allocate_fdof_to_data(n_fmdofs)
    # We use the conforming space dof ids to determine
    # the free dofs that constrain well-posed hanging 
    # dofs making them well posed (see Figure 3b of [R])
    _fill_fdof_to_is_agg!(fdof_to_is_agg,
                          acell_to_acellin,
                          acell_to_dof_ids) # On conforming space
    _fill_fdof_to_cell_and_ldof!(fdof_to_acell,
                                 fdof_to_ldof,
                                 acell_to_nc_dof_ids, # On non-conforming space
                                 acell_to_gcell,
                                 DOF_to_mDOFs)

    ##### END OLD

    ##### NEW

    # [!] fdofs follow the numbering of the non-conforming space
    dof_to_is_agg, dof_to_acell, dof_to_ldof = 
      _allocate_fdof_to_data(n_fdofs)
    # We use the conforming space dof ids to determine
    # the free dofs that constrain well-posed hanging 
    # dofs making them well posed (see Figure 3b of [R])

    # Arrays mapping dofs to hdofs already needed here
    hdof_to_dof = sDOF_to_dof
    dof_to_hdof = zeros(Int32,n_fdofs)
    dof_to_hdof[hdof_to_dof] .= 1:n_hdofs
    # This info can also be inferred when dof_to_hdof[dof] = 0
    dof_to_is_hdof = fill(false,n_fdofs)
    dof_to_is_hdof[hdof_to_dof] .= true
    hdof_to_dofs = sDOF_to_dofs

    _fill_dof_to_is_agg!(dof_to_is_agg,
                         acell_to_acellin,
                         acell_to_nc_dof_ids,
                         dof_to_is_hdof,
                         dof_to_hdof,
                         hdof_to_dofs) # On non-conforming space
    _fill_dof_to_cell_and_ldof!(dof_to_acell,
                                dof_to_ldof,
                                acell_to_nc_dof_ids, # On non-conforming space
                                acell_to_gcell)
    
    ##### END NEW

    # DEBUG
    dof_to_is_fdof = findall(.!dof_to_is_hdof)

    ##### OLD

    agg_fdof_to_fdof = findall(fdof_to_is_agg)
    fdof_to_agg_fdof = zeros(Int32,n_fmdofs)
    n_agg_fdofs = length(agg_fdof_to_fdof)
    fdof_to_agg_fdof[agg_fdof_to_fdof] .= 1:n_agg_fdofs

    # For refactor, need it before.
    #
    # hdof_to_dof = findall(length.(DOF_to_mDOFs).>1)
    # dof_to_hdof = zeros(Int32,n_fdofs)
    # dof_to_hdof[hdof_to_dof] .= 1:n_hdofs
    
    ##### END OLD

    ##### NEW

    agg_dof_to_dof = findall(dof_to_is_agg)
    dof_to_agg_dof = zeros(Int32,n_fdofs)
    n_agg_dofs = length(agg_dof_to_dof)
    dof_to_agg_dof[agg_dof_to_dof] .= 1:n_agg_dofs

    # Beware, different numbering of agg_dofs:
    # dof_to_agg_dof[dof_to_is_fdof] = Int32[1, 2, 3, 0, 4, 5, 0, 6, 7, 0, 0, 8, 9, 0, 10, 11, 12, 13, 0, 14, 15, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 0, 0, 0, 0, 34, 35, 36, 37, 38, 39, 40, 0, 0, 0, 41, 42, 43, 44, 0, 0, 45, 46, 47, 0, 0, 0, 60, 61, 62, 63, 64, 0, 65, 66, 67, 68, 0, 0, 69, 70, 0, 71, 72, 0, 0]
    # fdof_to_agg_fdof = Int32[1, 2, 3, 0, 4, 5, 0, 6, 7, 0, 0, 8, 9, 0, 10, 11, 12, 13, 0, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 0, 0, 0, 0, 28, 29, 30, 31, 32, 33, 34, 0, 0, 0, 35, 36, 37, 38, 0, 0, 39, 40, 41, 0, 0, 0, 42, 43, 44, 45, 46, 0, 47, 48, 49, 50, 0, 0, 51, 52, 0, 53, 54, 0, 0]

    ##### END NEW

    # [!] dofs follow the numbering of the non-conforming space
    agg_fdof_to_dof = collect(lazy_map(Reindex(mDOF_to_DOF),agg_fdof_to_fdof))
    aggdof_to_dof = vcat(agg_fdof_to_dof,hdof_to_dof)
    
    n_cdofs = n_agg_fdofs + n_hdofs # Total = ill-posed free + hanging dofs
    aggdof_to_dofs_ptrs = zeros(Int32,n_cdofs+1)

    _fill_aggdof_to_dofs_ptrs!(aggdof_to_dofs_ptrs,
                               agg_fdof_to_fdof, # On conforming space
                               fdof_to_acell,    # On conforming space
                               acell_to_acellin,
                               acell_to_dof_ids, # On conforming space
                               acell_to_is_owned)

    # A hanging dof might be mapped to >1 root cells
    hdof_to_is_agg = fill(true,n_hdofs)
    # We could further filter those that are 
    # only constrained by well-posed free dofs.
    # We do not do it right now, for simplicity.
    _fill_hdof_to_is_agg!(hdof_to_is_agg,
                          dof_to_hdof, # On non-conforming space
                          acell_to_acellin,
                          acell_to_nc_dof_ids,
                          DOF_to_mDOFs)
    @show hdof_to_is_agg

    # In distributed, how can I determine the unique number  
    # of constraining dofs per ill-posed hanging dof?
    _fill_hdof_to_dofs_ptrs!(aggdof_to_dofs_ptrs,
                             n_fdofs,
                             n_fmdofs,
                             n_agg_fdofs,
                             fdof_to_is_agg,
                             hdof_to_is_agg,
                             hdof_to_dof,
                             fdof_to_agg_fdof,
                             DOF_to_mDOFs)

    aggdof_to_dofs_data, aggdof_to_coeffs_data = 
      _allocate_aggdof_to_data(aggdof_to_dofs_ptrs,
                               acell_to_coeffs)

    _fill_aggdof_to_dofs_data!(aggdof_to_dofs_data,
                               aggdof_to_dofs_ptrs,
                               agg_fdof_to_fdof,
                               fdof_to_acell,
                               acell_to_acellin,
                               acell_to_dof_ids, # On conforming space
                               acell_to_is_owned,
                               mDOF_to_DOF)

    _fill_aggdof_to_coeffs_data!(aggdof_to_coeffs_data,
                                 aggdof_to_dofs_ptrs,
                                 agg_fdof_to_fdof,
                                 fdof_to_acell,
                                 fdof_to_ldof,
                                 acellin_constraints,
                                 acell_to_coeffs,
                                 acell_to_proj,
                                 acell_to_is_owned)

    _fill_hdof_to_dofs_data!(aggdof_to_dofs_data,
                             aggdof_to_dofs_ptrs,
                             n_fdofs,
                             n_fmdofs,
                             n_agg_fdofs,
                             fdof_to_is_agg,
                             hdof_to_is_agg,
                             hdof_to_dof,
                             fdof_to_agg_fdof,
                             DOF_to_mDOFs,
                             mDOF_to_DOF)

    _fill_hdof_to_coeffs_data!(aggdof_to_coeffs_data,
                               aggdof_to_dofs_ptrs,
                               n_fmdofs,
                               n_agg_fdofs,
                               fdof_to_is_agg,
                               hdof_to_is_agg,
                               hdof_to_dof,
                               fdof_to_agg_fdof,
                               DOF_to_mDOFs,
                               DOF_to_coeffs)

    aggdof_to_dofs   = Table(aggdof_to_dofs_data,  
                             aggdof_to_dofs_ptrs)
    aggdof_to_coeffs = Table(aggdof_to_coeffs_data,
                             aggdof_to_dofs_ptrs)

    aggdof_to_dof, aggdof_to_dofs, aggdof_to_coeffs
  end

  function _fill_fdof_to_is_agg!(
    fdof_to_is_agg,
    acell_to_acellin,
    acell_to_dof_ids)

    cache = array_cache(acell_to_dof_ids)
    for (acell,acellin) in enumerate(acell_to_acellin)
      iscut = acell != acellin
      if !iscut
        dofs = getindex!(cache,acell_to_dof_ids,acell)
        for dof in dofs
          (dof > 0) && (fdof_to_is_agg[dof] = false)
        end
      end
    end
    nothing
  end

  function _fill_dof_to_is_agg!(
    dof_to_is_agg,
    acell_to_acellin,
    acell_to_dof_ids,
    dof_to_is_hdof,
    dof_to_hdof,
    hdof_to_dofs)

    cache = array_cache(acell_to_dof_ids)
    for (acell,acellin) in enumerate(acell_to_acellin)
      iscut = acell != acellin
      if !iscut
        dofs = getindex!(cache,acell_to_dof_ids,acell)
        for dof in dofs
          dof < 0 && continue
          dof_to_is_agg[dof] = false
          if dof_to_is_hdof[dof]
            # Mark as well-posed all constraining free dofs
            for mdof in hdof_to_dofs[dof_to_hdof[dof]]
              (mdof > 0) && (dof_to_is_agg[mdof] = false)
            end
          end
        end
      end
    end
    nothing
  end

  function _fill_fdof_to_cell_and_ldof!(
    fdof_to_acell,
    fdof_to_ldof,
    acell_to_nc_dof_ids,
    acell_to_gcell,
    DOF_to_mDOFs)

    cache = array_cache(acell_to_nc_dof_ids)
    for (acell,gcell) in enumerate(acell_to_gcell)
      dofs = getindex!(cache,acell_to_nc_dof_ids,acell)
      for (ldof,dof) in enumerate(dofs)
        if (dof > 0) && (length(DOF_to_mDOFs[dof]) == 1) # It's a free dof
          fdof = DOF_to_mDOFs[dof][1]
          acell_dof = fdof_to_acell[fdof]
          if acell_dof == 0 || (gcell > acell_to_gcell[acell_dof])
            fdof_to_acell[fdof] = acell
            fdof_to_ldof[fdof] = ldof
          end
        end
      end
    end
    nothing
  end

  function _fill_dof_to_cell_and_ldof!(
    dof_to_acell,
    dof_to_ldof,
    acell_to_nc_dof_ids,
    acell_to_gcell)

    cache = array_cache(acell_to_nc_dof_ids)
    for (acell,gcell) in enumerate(acell_to_gcell)
      dofs = getindex!(cache,acell_to_nc_dof_ids,acell)
      for (ldof,dof) in enumerate(dofs)
        if (dof > 0)
          acell_dof = dof_to_acell[dof]
          if acell_dof == 0 || (gcell > acell_to_gcell[acell_dof])
            dof_to_acell[dof] = acell
            dof_to_ldof[dof] = ldof
          end
        end
      end
    end
    nothing
  end

  function _allocate_hdof_to_data(n_hdofs)
    hdof_to_is_agg = fill(true,n_hdofs)
    hdof_to_dof = zeros(Int32,n_hdofs)
    hdof_to_is_agg, hdof_to_dof
  end

  function _fill_hdof_to_is_agg!(hdof_to_is_agg,
                                 dof_to_hdof, # On non-conforming space
                                 acell_to_acellin,
                                 acell_to_nc_dof_ids,
                                 DOF_to_mDOFs)
    
    cache = array_cache(acell_to_nc_dof_ids)
    for (acell,acellin) in enumerate(acell_to_acellin)
      iscut = acell != acellin
      if !iscut
        dofs = getindex!(cache,acell_to_nc_dof_ids,acell)
        for dof in dofs
          if (dof > 0) && (length(DOF_to_mDOFs[dof])>1) # It's a hanging dof
            hdof_to_is_agg[dof_to_hdof[dof]] = false
          end
        end
      end
    end
    nothing
  end

  function _fill_hdof_to_dofs_ptrs!(
    aggdof_to_dofs_ptrs,
    n_fdofs,
    n_fmdofs,
    n_agg_fdofs,
    fdof_to_is_agg,
    hdof_to_is_agg,
    hdof_to_dof,
    fdof_to_agg_fdof,
    DOF_to_mDOFs)

    for (hdof,hdof_is_agg) in enumerate(hdof_to_is_agg)
      if hdof_is_agg # ill-posed
        DOF = hdof_to_dof[hdof]
        for mDOF in DOF_to_mDOFs[DOF]
          if mDOF <= n_fmdofs # free master dof
            fdof = mDOF
            fdof_is_agg = fdof_to_is_agg[fdof]
            if fdof_is_agg # ill-posed
              f_aggdof = fdof_to_agg_fdof[fdof]
              ndofs_f_aggdof = aggdof_to_dofs_ptrs[f_aggdof+1]
              aggdof_to_dofs_ptrs[n_agg_fdofs+hdof+1] += ndofs_f_aggdof
            else # well-posed
              aggdof_to_dofs_ptrs[n_agg_fdofs+hdof+1] += 1
            end 
          else # Dirichlet master dof
            aggdof_to_dofs_ptrs[n_agg_fdofs+hdof+1] += 1
          end
        end
      else # well-posed
        aggdof_to_dofs_ptrs[n_agg_fdofs+hdof+1] = 
          length(DOF_to_mDOFs[hdof_to_dof[hdof]])
      end
    end
    nothing
  end

  function _fill_aggdof_to_dofs_data!(
    aggdof_to_dofs_data,
    aggdof_to_dofs_ptrs,
    agg_fdof_to_fdof,
    fdof_to_acell,
    acell_to_acellin,
    acell_to_dof_ids,
    acell_to_is_owned,
    mDOF_to_DOF)

    cache = array_cache(acell_to_dof_ids)
    for (aggdof,fdof) in enumerate(agg_fdof_to_fdof)
      acell = fdof_to_acell[fdof]
      ! acell_to_is_owned[acell] && continue
      acellin = acell_to_acellin[acell]
      dofs = getindex!(cache,acell_to_dof_ids,acellin)
      p = aggdof_to_dofs_ptrs[aggdof]-1
      for (i,dof) in enumerate(dofs)
        dof = dof > 0 ? mDOF_to_DOF[dof] : dof
        aggdof_to_dofs_data[p+i] = dof
      end
    end
    nothing
  end

  function _fill_aggdof_to_coeffs_data!(
    aggdof_to_coeffs_data,
    aggdof_to_dofs_ptrs,
    agg_fdof_to_fdof,
    fdof_to_acell,
    fdof_to_ldof,
    acellin_constraints,
    acell_to_coeffs,
    acell_to_proj,
    acell_to_is_owned)

    cache2 = array_cache(acell_to_coeffs)
    cache3 = array_cache(acell_to_proj)
    cache4 = array_cache(acellin_constraints)

    T = eltype(eltype(acell_to_coeffs))
    z = zero(T)

    for (aggdof,fdof) in enumerate(agg_fdof_to_fdof)
      acell = fdof_to_acell[fdof]
      ! acell_to_is_owned[acell] && continue
      coeffs = getindex!(cache2,acell_to_coeffs,acell)
      proj = getindex!(cache3,acell_to_proj,acell)
      constr = getindex!(cache4,acellin_constraints,acell) # lmdof x ldof
      ldof = fdof_to_ldof[fdof]
      p = aggdof_to_dofs_ptrs[aggdof]-1
      for lmdof in 1:size(constr,1) # lmdof on acellin
        coeff = z
        for b in 1:size(proj,2)
          for c in 1:size(coeffs,2)
            coeff += coeffs[ldof,c]*proj[c,b]*constr[lmdof,b]
          end
        end
        aggdof_to_coeffs_data[p+lmdof] = coeff
      end
    end
    nothing
  end

  function _fill_hdof_to_dofs_data!(
    aggdof_to_dofs_data,
    aggdof_to_dofs_ptrs,
    n_fdofs,
    n_fmdofs,
    n_agg_fdofs,
    fdof_to_is_agg,
    hdof_to_is_agg,
    hdof_to_dof,
    fdof_to_agg_fdof,
    DOF_to_mDOFs,
    mDOF_to_DOF)

    for (hdof,hdof_is_agg) in enumerate(hdof_to_is_agg)
      if hdof_is_agg # ill-posed
        dofs = eltype(aggdof_to_dofs_data)[]
        DOF = hdof_to_dof[hdof]
        for mDOF in DOF_to_mDOFs[DOF]
          if mDOF <= n_fmdofs # free master dof
            fdof = mDOF
            fdof_is_agg = fdof_to_is_agg[fdof]
            if fdof_is_agg # ill-posed
              f_aggdof = fdof_to_agg_fdof[fdof]
              _ini = aggdof_to_dofs_ptrs[f_aggdof]
              _end = aggdof_to_dofs_ptrs[f_aggdof+1]-1
              _dofs = aggdof_to_dofs_data[_ini:_end]
              dofs = vcat(dofs,_dofs)
            else # well-posed
              dofs = vcat(dofs,[mDOF_to_DOF[fdof]])
            end 
          else # Dirichlet master dof
            ddof = -(mDOF_to_DOF[mDOF]-n_fdofs)
            dofs = vcat(dofs,[ddof])
          end
        end
      else # well-posed
        dofs = mDOF_to_DOF[DOF_to_mDOFs[hdof_to_dof[hdof]]]
      end
      p = aggdof_to_dofs_ptrs[n_agg_fdofs+hdof]-1
      for (i,dof) in enumerate(dofs)
        aggdof_to_dofs_data[p+i] = dof
      end
    end
    nothing
  end

  function _fill_hdof_to_coeffs_data!(
    aggdof_to_coeffs_data,
    aggdof_to_dofs_ptrs,
    n_fmdofs,
    n_agg_fdofs,
    fdof_to_is_agg,
    hdof_to_is_agg,
    hdof_to_dof,
    fdof_to_agg_fdof,
    DOF_to_mDOFs,
    DOF_to_coeffs)

    for (hdof,hdof_is_agg) in enumerate(hdof_to_is_agg)
      if hdof_is_agg # ill-posed
        coeffs = eltype(aggdof_to_coeffs_data)[]
        DOF = hdof_to_dof[hdof]
        for (i,mDOF) in enumerate(DOF_to_mDOFs[DOF])
          if mDOF <= n_fmdofs # free master dof
            fdof = mDOF
            fdof_is_agg = fdof_to_is_agg[fdof]
            if fdof_is_agg # ill-posed
              f_aggdof = fdof_to_agg_fdof[fdof]
              _ini = aggdof_to_dofs_ptrs[f_aggdof]
              _end = aggdof_to_dofs_ptrs[f_aggdof+1]-1
              _coeffs = DOF_to_coeffs[DOF][i] .* 
                aggdof_to_coeffs_data[_ini:_end]
              coeffs = vcat(coeffs,_coeffs)
            else # well-posed
              coeffs = vcat(coeffs,[DOF_to_coeffs[DOF][i]])
            end
          else # Dirichlet master dof
            coeffs = vcat(coeffs,[DOF_to_coeffs[DOF][i]])
          end
        end
      else # well-posed
        DOF = hdof_to_dof[hdof]
        coeffs = DOF_to_coeffs[DOF]
      end
      p = aggdof_to_dofs_ptrs[n_agg_fdofs+hdof]-1
      for (i,coeff) in enumerate(coeffs)
        aggdof_to_coeffs_data[p+i] = coeff
      end
    end
    nothing
  end

end