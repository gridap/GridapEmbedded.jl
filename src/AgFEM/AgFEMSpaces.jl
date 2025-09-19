
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
  acell_to_gcell,
  acell_to_is_owned=fill(true,length(acell_to_acellin)))

  fdof_to_is_agg, fdof_to_acell, fdof_to_ldof = 
    _allocate_fdof_to_data(n_fdofs)

  _fill_fdof_to_data!(fdof_to_is_agg,
                      fdof_to_acell,
                      fdof_to_ldof,
                      acell_to_acellin,
                      acell_to_dof_ids,
                      acell_to_gcell)

  aggdof_to_fdof = findall(fdof_to_is_agg)

  n_aggdofs = length(aggdof_to_fdof)
  aggdof_to_dofs_ptrs = zeros(Int32,n_aggdofs+1)

  _fill_aggdof_to_dofs_ptrs!(aggdof_to_dofs_ptrs,
                             aggdof_to_fdof,
                             fdof_to_acell,
                             acell_to_acellin,
                             acell_to_dof_ids,
                             acell_to_is_owned)

  aggdof_to_dofs_data, aggdof_to_coeffs_data = 
    _allocate_aggdof_to_data(aggdof_to_dofs_ptrs,
                             acell_to_coeffs)

  _fill_aggdof_to_dofs_data!(aggdof_to_dofs_data,
                             aggdof_to_dofs_ptrs,
                             aggdof_to_fdof,
                             fdof_to_acell,
                             acell_to_acellin,
                             acell_to_dof_ids,
                             acell_to_is_owned)

  _fill_aggdof_to_coeffs_data!(aggdof_to_coeffs_data,
                               aggdof_to_dofs_ptrs,
                               aggdof_to_fdof,
                               fdof_to_acell,
                               fdof_to_ldof,
                               acell_to_coeffs,
                               acell_to_proj,
                               acell_to_is_owned)

  aggdof_to_dofs   = Table(aggdof_to_dofs_data,  
                           aggdof_to_dofs_ptrs)
  aggdof_to_coeffs = Table(aggdof_to_coeffs_data,
                           aggdof_to_dofs_ptrs)

  aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs
end

function _allocate_fdof_to_data(n_fdofs)
  fdof_to_is_agg = fill(true,n_fdofs)
  fdof_to_acell = zeros(Int32,n_fdofs)
  fdof_to_ldof = zeros(Int16,n_fdofs)
  fdof_to_is_agg, fdof_to_acell, fdof_to_ldof
end

function _fill_fdof_to_data!(
  fdof_to_is_agg,
  fdof_to_acell,
  fdof_to_ldof,
  acell_to_acellin,
  acell_to_dof_ids,
  acell_to_gcell)

  # RMK: There can be owned dofs sitting around non-owned 
  # cells that are mapped to an acell for which acellin is
  # zero (i.e., outside the ghost layer).
  cache = array_cache(acell_to_dof_ids)
  for (acell, acellin) in enumerate(acell_to_acellin)
    iscut = acell != acellin
    dofs = getindex!(cache,acell_to_dof_ids,acell)
    gcell = acell_to_gcell[acell]
    for (ldof,dof) in enumerate(dofs)
      if dof > 0
        fdof = dof
        acell_dof = fdof_to_acell[fdof]
        fdof_to_is_agg[fdof] &= iscut
        if acell_dof == 0 || gcell > acell_to_gcell[acell_dof]
          fdof_to_acell[fdof] = acell
          fdof_to_ldof[fdof] = ldof
         end
      end
    end
  end
  nothing
end

function _fill_aggdof_to_dofs_ptrs!(
  aggdof_to_dofs_ptrs,
  aggdof_to_fdof,
  fdof_to_acell,
  acell_to_acellin,
  acell_to_dof_ids,
  acell_to_is_owned)

  cache = array_cache(acell_to_dof_ids)
  for (aggdof,fdof) in enumerate(aggdof_to_fdof)
    acell = fdof_to_acell[fdof]
    ! acell_to_is_owned[acell] && continue
    acellin = acell_to_acellin[acell]
    dofs = getindex!(cache,acell_to_dof_ids,acellin)
    aggdof_to_dofs_ptrs[aggdof+1] = length(dofs)
  end
  nothing
end 

function _allocate_aggdof_to_data(aggdof_to_dofs_ptrs,
                                  acell_to_coeffs)
  length_to_ptrs!(aggdof_to_dofs_ptrs)
  ndata = aggdof_to_dofs_ptrs[end]-1
  aggdof_to_dofs_data = zeros(Int,ndata)
  T = eltype(eltype(acell_to_coeffs))
  aggdof_to_coeffs_data = zeros(T,ndata)
  aggdof_to_dofs_data, aggdof_to_coeffs_data
end

function _fill_aggdof_to_dofs_data!(
  aggdof_to_dofs_data,
  aggdof_to_dofs_ptrs,
  aggdof_to_fdof,
  fdof_to_acell,
  acell_to_acellin,
  acell_to_dof_ids,
  acell_to_is_owned)

  cache = array_cache(acell_to_dof_ids)
  for (aggdof,fdof) in enumerate(aggdof_to_fdof)
    acell = fdof_to_acell[fdof]
    ! acell_to_is_owned[acell] && continue
    acellin = acell_to_acellin[acell]
    dofs = getindex!(cache,acell_to_dof_ids,acellin)
    p = aggdof_to_dofs_ptrs[aggdof]-1
    for (i,dof) in enumerate(dofs)
      aggdof_to_dofs_data[p+i] = dof
    end
  end
  nothing
end

function _fill_aggdof_to_coeffs_data!(
  aggdof_to_coeffs_data,
  aggdof_to_dofs_ptrs,
  aggdof_to_fdof,
  fdof_to_acell,
  fdof_to_ldof,
  acell_to_coeffs,
  acell_to_proj,
  acell_to_is_owned)

  cache2 = array_cache(acell_to_coeffs)
  cache3 = array_cache(acell_to_proj)

  T = eltype(eltype(acell_to_coeffs))
  z = zero(T)

  for (aggdof,fdof) in enumerate(aggdof_to_fdof)
    acell = fdof_to_acell[fdof]
    ! acell_to_is_owned[acell] && continue
    coeffs = getindex!(cache2,acell_to_coeffs,acell)
    proj = getindex!(cache3,acell_to_proj,acell)
    ldof = fdof_to_ldof[fdof]
    p = aggdof_to_dofs_ptrs[aggdof]-1
    for b in 1:size(proj,2)
      coeff = z
      for c in 1:size(coeffs,2)
        coeff += coeffs[ldof,c]*proj[c,b]
      end
      aggdof_to_coeffs_data[p+b] = coeff
    end
  end
  nothing
end

###########################################################################################
# Unused functions

function _alloc_and_fill_aggdof_to_dofs_ptrs(
  n_fdofs,
  acell_to_acellin,
  acell_to_dof_ids,
  acell_to_gcell,
  acell_to_is_owned)

  fdof_to_is_agg, fdof_to_acell, fdof_to_ldof = 
    _allocate_fdof_to_data(n_fdofs)

  _fill_fdof_to_data!(fdof_to_is_agg,
                      fdof_to_acell,
                      fdof_to_ldof,
                      acell_to_acellin,
                      acell_to_dof_ids,
                      acell_to_gcell)

  aggdof_to_fdof = findall(fdof_to_is_agg)

  n_aggdofs = length(aggdof_to_fdof)
  aggdof_to_dofs_ptrs = zeros(Int32,n_aggdofs+1)

  _fill_aggdof_to_dofs_ptrs!(aggdof_to_dofs_ptrs,
                             aggdof_to_fdof,
                             fdof_to_acell,
                             acell_to_acellin,
                             acell_to_dof_ids,
                             acell_to_is_owned)

  aggdof_to_dofs_ptrs, aggdof_to_fdof, fdof_to_acell, fdof_to_ldof
end

function _alloc_and_fill_aggdof_to_dofs_data(
  aggdof_to_fdof,
  aggdof_to_dofs_ptrs,
  acell_to_acellin,
  acell_to_dof_ids,
  acell_to_coeffs,
  acell_to_proj,
  acell_to_is_owned,
  fdof_to_acell,
  fdof_to_ldof)

  aggdof_to_dofs_data, aggdof_to_coeffs_data = 
    _allocate_aggdof_to_data(aggdof_to_dofs_ptrs,
                             acell_to_coeffs)

  _fill_aggdof_to_dofs_data!(aggdof_to_dofs_data,
                             aggdof_to_dofs_ptrs,
                             aggdof_to_fdof,
                             fdof_to_acell,
                             acell_to_acellin,
                             acell_to_dof_ids,
                             acell_to_is_owned)

  _fill_aggdof_to_coeffs_data!(aggdof_to_coeffs_data,
                               aggdof_to_dofs_ptrs,
                               aggdof_to_fdof,
                               fdof_to_acell,
                               fdof_to_ldof,
                               acell_to_coeffs,
                               acell_to_proj,
                               acell_to_is_owned)

  aggdof_to_dofs   = Table(aggdof_to_dofs_data,  
                           aggdof_to_dofs_ptrs)
  aggdof_to_coeffs = Table(aggdof_to_coeffs_data,
                           aggdof_to_dofs_ptrs)

  aggdof_to_dofs, aggdof_to_coeffs
end

function _fill_acell_to_acellin_and_to_gcell(trian_a,bgcell_to_bgcellin,bgcell_to_gcell)
  D = num_cell_dims(trian_a)
  glue = get_glue(trian_a,Val(D))
  acell_to_bgcell = glue.tface_to_mface
  bgcell_to_acell = glue.mface_to_tface
  acell_to_bgcellin = 
    collect(lazy_map(Reindex(bgcell_to_bgcellin),acell_to_bgcell))
  T = eltype(bgcell_to_acell)
  acell_to_acellin = map(acell_to_bgcellin) do bgcin
    iszero(bgcin) ? zero(T) : T(bgcell_to_acell[bgcin])
  end
  acell_to_gcell = lazy_map(Reindex(bgcell_to_gcell),acell_to_bgcell)
  acell_to_acellin, acell_to_gcell
end
