
function AgFEMSpace(
  f::SingleFieldFESpace,
  bgcell_to_bgcellin::AbstractVector,
  g::SingleFieldFESpace=f,
  args...)

  @assert get_triangulation(f) === get_triangulation(g)
  AgFEMSpace(f,bgcell_to_bgcellin,get_fe_basis(g),get_fe_dof_basis(g),args...)
end

# Note: cell is in fact bgcell in this function since f will usually be an ExtendedFESpace
function AgFEMSpace(
  f::SingleFieldFESpace,
  bgcell_to_bgcellin::AbstractVector,
  shfns_g::CellField,
  dofs_g::CellDof,
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

  aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs = _setup_agfem_constraints(
    num_free_dofs(f),
    acell_to_acellin,
    acell_to_dof_ids,
    acell_to_coeffs,
    acell_to_proj,
    acell_to_gcell)

  FESpaceWithLinearConstraints(aggdof_to_fdof,aggdof_to_dofs,aggdof_to_coeffs,f)
end

function _setup_agfem_constraints(
  n_fdofs,
  acell_to_acellin,
  acell_to_dof_ids,
  acell_to_coeffs,
  acell_to_proj,
  acell_to_gcell)

  n_acells = length(acell_to_acellin)
  fdof_to_isagg = fill(true,n_fdofs)
  fdof_to_acell = zeros(Int32,n_fdofs)
  fdof_to_ldof = zeros(Int16,n_fdofs)
  cache = array_cache(acell_to_dof_ids)
  for acell in 1:n_acells
    acellin = acell_to_acellin[acell]
    iscut = acell != acellin
    dofs = getindex!(cache,acell_to_dof_ids,acell)
    gcell = acell_to_gcell[acell]
    for (ldof,dof) in enumerate(dofs)
      if dof > 0
        fdof = dof
        acell_dof = fdof_to_acell[fdof]
        if acell_dof == 0 || gcell > acell_to_gcell[acell_dof]
          fdof_to_acell[fdof] = acell
          fdof_to_isagg[fdof] = iscut && fdof_to_isagg[fdof]
          fdof_to_ldof[fdof] = ldof
         end
      end
    end
  end

  aggdof_to_fdof = findall(fdof_to_isagg)

  n_aggdofs = length(aggdof_to_fdof)
  aggdof_to_dofs_ptrs = zeros(Int32,n_aggdofs+1)

  for aggdof in 1:n_aggdofs
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    acellin = acell_to_acellin[acell]
    dofs = getindex!(cache,acell_to_dof_ids,acellin)
    aggdof_to_dofs_ptrs[aggdof+1] = length(dofs)
  end

  length_to_ptrs!(aggdof_to_dofs_ptrs)
  ndata = aggdof_to_dofs_ptrs[end]-1
  aggdof_to_dofs_data = zeros(Int,ndata)

  for aggdof in 1:n_aggdofs
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    acellin = acell_to_acellin[acell]
    dofs = getindex!(cache,acell_to_dof_ids,acellin)
    p = aggdof_to_dofs_ptrs[aggdof]-1
    for (i,dof) in enumerate(dofs)
      aggdof_to_dofs_data[p+i] = dof
    end
  end

  aggdof_to_dofs = Table(aggdof_to_dofs_data,aggdof_to_dofs_ptrs)

  cache2 = array_cache(acell_to_coeffs)
  cache3 = array_cache(acell_to_proj)

  T = eltype(eltype(acell_to_coeffs))
  z = zero(T)

  aggdof_to_coefs_data = zeros(T,ndata)
  for aggdof in 1:n_aggdofs
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    coeffs = getindex!(cache2,acell_to_coeffs,acell)
    proj = getindex!(cache3,acell_to_proj,acell)
    ldof = fdof_to_ldof[fdof]
    p = aggdof_to_dofs_ptrs[aggdof]-1
    for b in 1:size(proj,2)
      coeff = z
      for c in 1:size(coeffs,2)
        coeff += coeffs[ldof,c]*proj[c,b]
      end
      aggdof_to_coefs_data[p+b] = coeff
    end
  end

  aggdof_to_coeffs = Table(aggdof_to_coefs_data,aggdof_to_dofs_ptrs)

  aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs
end
