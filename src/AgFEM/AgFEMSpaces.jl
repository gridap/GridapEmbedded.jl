
function AgFEMSpace(f::ExtendedFESpace,acell_to_cellin::AbstractVector)

  trian = f.trian
  acell_to_cell = trian.cell_to_oldcell
  cell_to_acell = trian.oldcell_to_cell
  acell_to_dofs = reindex(get_cell_dofs(f),acell_to_cell)
  n_fdofs = num_free_dofs(f) 
  acell_to_basis = reindex(get_cell_basis(f),acell_to_cellin)
  acell_to_dof_basis = reindex(get_cell_dof_basis(f),acell_to_cell)
  @notimplementedif is_in_ref_space(acell_to_basis)
  @notimplementedif is_in_ref_space(acell_to_dof_basis)
  acell_to_coeffs = evaluate(acell_to_dof_basis,acell_to_basis)

  aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs = _setup_agfem_constraints(
    n_fdofs,
    acell_to_cellin,
    acell_to_cell,
    cell_to_acell,
    acell_to_dofs,
    acell_to_coeffs)

  FESpaceWithLinearConstraints(aggdof_to_fdof,aggdof_to_dofs,aggdof_to_coeffs,f)
end

function _setup_agfem_constraints(
  n_fdofs,
  acell_to_cellin,
  acell_to_cell,
  cell_to_acell,
  acell_to_dofs,
  acell_to_coeffs)

  n_acells = length(acell_to_cell)
  fdof_to_isagg = fill(true,n_fdofs)
  fdof_to_acell = zeros(Int32,n_fdofs)
  fdof_to_ldof = zeros(Int8,n_fdofs)
  cache = array_cache(acell_to_dofs)
  for acell in 1:n_acells
    iscut = acell_to_cell[acell] != acell_to_cellin[acell]
    dofs = getindex!(cache,acell_to_dofs,acell)
    cellin = acell_to_cellin[acell]
    for (ldof,dof) in enumerate(dofs)
      if dof > 0
        fdof = dof
        fdof_to_isagg[fdof] = iscut && fdof_to_isagg[fdof]
        fdof_to_acell[fdof] = acell
        fdof_to_ldof[fdof] = ldof
      end
    end
  end

  aggdof_to_fdof = findall(fdof_to_isagg)

  n_aggdofs = length(aggdof_to_fdof)
  aggdof_to_dofs_ptrs = zeros(Int32,n_aggdofs+1)

  for aggdof in 1:n_aggdofs
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    cellin = acell_to_cellin[acell]
    acell = cell_to_acell[cellin]
    dofs = getindex!(cache,acell_to_dofs,acell)
    aggdof_to_dofs_ptrs[aggdof+1] = length(dofs)
  end

  length_to_ptrs!(aggdof_to_dofs_ptrs)
  ndata = aggdof_to_dofs_ptrs[end]-1
  aggdof_to_dofs_data = zeros(Int,ndata)

  for aggdof in 1:n_aggdofs
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    cellin = acell_to_cellin[acell]
    acell = cell_to_acell[cellin]
    dofs = getindex!(cache,acell_to_dofs,acell)
    p = aggdof_to_dofs_ptrs[aggdof]-1
    for (i,dof) in enumerate(dofs)
      aggdof_to_dofs_data[p+i] = dof
    end
  end

  aggdof_to_dofs = Table(aggdof_to_dofs_data,aggdof_to_dofs_ptrs)

  cache2 = array_cache(acell_to_coeffs)

  aggdof_to_coefs_data = zeros(Float64,ndata)
  for aggdof in 1:n_aggdofs
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    coeffs = getindex!(cache2,acell_to_coeffs,acell)
    ldof = fdof_to_ldof[fdof]
    p = aggdof_to_dofs_ptrs[aggdof]-1
    n_coeffs = size(coeffs,2)
    for i in 1:n_coeffs
      coeff = coeffs[ldof,i]
      aggdof_to_coefs_data[p+i] = coeff
    end
  end

  aggdof_to_coeffs = Table(aggdof_to_coefs_data,aggdof_to_dofs_ptrs)

  aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs
end

