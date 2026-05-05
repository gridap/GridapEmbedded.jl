
using Gridap
using PartitionedArrays
using GridapDistributed

using Gridap.Arrays
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData
using Gridap.Fields

using Gridap.Geometry: get_active_model

using GridapEmbedded

using BenchmarkTools
using FillArrays
using LinearAlgebra

############################################################################################
# Cartesian map optimization
# Note: Not really that big of a difference. 

struct InverseCartesianMap{D,T,L} <: AbstractArray{AffineField{D,D,T,L},D}
  data::CartesianDescriptor{D,T,typeof(identity)}
  function InverseCartesianMap(des::CartesianDescriptor{D,T}) where {D,T}
    L = D*D
    new{D,T,L}(des)
  end
end

Base.size(a::InverseCartesianMap) = a.data.partition

Base.IndexStyle(::Type{<:InverseCartesianMap}) = IndexCartesian()

function Base.getindex(a::InverseCartesianMap{D,T},I::Vararg{Integer,D}) where {D,T}
  x0 = a.data.origin
  dx = a.data.sizes
  inv_dx = map(inv, dx)
  fwd_origin = x0 + Point(ntuple(d -> T(I[d]-1)*dx[d], D))
  inv_grad = diagonal_tensor(VectorValue(inv_dx))
  inv_origin = Point(ntuple(d -> -fwd_origin[d]*inv_dx[d], D))
  return AffineField(inv_grad, inv_origin)
end

function Arrays.lazy_map(k::typeof(Fields.inverse_map),a::Geometry.CartesianMap)
  InverseCartesianMap(a.data)
end

function Arrays.lazy_map(k::typeof(Fields.inverse_map),a::LazyArray{<:Fill{<:Reindex{<:Geometry.CartesianMap}}})
  println("Flag 0")
  imaps = lazy_map(inverse_map, a.maps.value.values)
  lazy_map(Reindex(imaps), a.args[1])
end

#######################

function setup_aggdata_v0(cutgeo,f,g)
  bgcell_to_bgcellin = aggregate(AggregateCutCellsByThreshold(0.8),cutgeo)

  shfns_g = get_fe_basis(g)
  dofs_g = get_fe_dof_basis(g)
  bgcell_to_gcell=1:length(bgcell_to_bgcellin)
  trian_a = get_triangulation(f)

  D = num_cell_dims(trian_a)
  glue = get_glue(trian_a,Val(D))
  acell_to_bgcell = glue.tface_to_mface
  bgcell_to_acell = glue.mface_to_tface
  acell_to_bgcellin = collect(lazy_map(Reindex(bgcell_to_bgcellin),acell_to_bgcell))
  acell_to_acellin = collect(lazy_map(Reindex(bgcell_to_acell),acell_to_bgcellin))
  acell_to_gcell = lazy_map(Reindex(bgcell_to_gcell),acell_to_bgcell)

  acell_phys_shapefuns_g = get_array(change_domain(shfns_g,PhysicalDomain()))
  acell_phys_root_shapefuns_g = lazy_map(Reindex(acell_phys_shapefuns_g),acell_to_acellin)
  root_shfns_g = GenericCellField(acell_phys_root_shapefuns_g,trian_a,PhysicalDomain())

  n_fdofs = num_free_dofs(f)
  dofs_f = get_fe_dof_basis(f)
  shfns_f = get_fe_basis(f)
  acell_to_coeffs = dofs_f(root_shfns_g)
  acell_to_proj = dofs_g(shfns_f)
  acell_to_dof_ids = get_cell_dof_ids(f)

  acell_to_is_owned=fill(true,length(acell_to_acellin))
  fdof_to_is_agg, fdof_to_acell, fdof_to_ldof = GridapEmbedded.AgFEM._allocate_fdof_to_data(n_fdofs)
  GridapEmbedded.AgFEM._fill_fdof_to_data!(
    fdof_to_is_agg,fdof_to_acell,fdof_to_ldof,acell_to_acellin,acell_to_dof_ids,acell_to_gcell
  )
  aggdof_to_fdof = findall(fdof_to_is_agg)
  n_aggdofs = length(aggdof_to_fdof)
  aggdof_to_dofs_ptrs = zeros(Int32,n_aggdofs+1)
  GridapEmbedded.AgFEM._fill_aggdof_to_dofs_ptrs!(
    aggdof_to_dofs_ptrs,aggdof_to_fdof,fdof_to_acell,acell_to_acellin,
    acell_to_dof_ids,acell_to_is_owned
  )
  aggdof_to_dofs_data, aggdof_to_coeffs_data = GridapEmbedded.AgFEM._allocate_aggdof_to_data(
      aggdof_to_dofs_ptrs,acell_to_coeffs
  )
  GridapEmbedded.AgFEM._fill_aggdof_to_dofs_data!(
    aggdof_to_dofs_data,aggdof_to_dofs_ptrs,aggdof_to_fdof,fdof_to_acell,
    acell_to_acellin,acell_to_dof_ids,acell_to_is_owned
  )
  GridapEmbedded.AgFEM._fill_aggdof_to_coeffs_data!(
    aggdof_to_coeffs_data, aggdof_to_dofs_ptrs, aggdof_to_fdof, fdof_to_acell, fdof_to_ldof, 
    acell_to_coeffs, acell_to_proj, acell_to_is_owned
  )
  aggdof_to_dofs   = Table(aggdof_to_dofs_data,  aggdof_to_dofs_ptrs)
  aggdof_to_coeffs = Table(aggdof_to_coeffs_data,aggdof_to_dofs_ptrs)

  return (; 
    aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs, fdof_to_acell, 
    fdof_to_ldof, acell_to_gcell, acell_to_coeffs, acell_to_proj, acell_to_is_owned
  )
end

function _fill_aggdof_to_coeffs_data_v0!(
  aggdof_to_coeffs_data,
  aggdof_to_dofs_ptrs,
  aggdof_to_fdof,
  fdof_to_acell,
  fdof_to_ldof,
  acell_to_coeffs,
  acell_to_proj,
  acell_to_is_owned,
  perm = eachindex(aggdof_to_fdof)
)
  cache2 = array_cache(acell_to_coeffs)
  cache3 = array_cache(acell_to_proj)

  T = eltype(eltype(acell_to_coeffs))
  z = zero(T)

  for aggdof in perm
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    ! acell_to_is_owned[acell] && continue
    coeffs = getindex!(cache2,acell_to_coeffs,acell)
    proj = getindex!(cache3,acell_to_proj,acell)
    ldof = fdof_to_ldof[fdof]
    p = aggdof_to_dofs_ptrs[aggdof]-1
    for b in axes(proj,2)
      coeff = z
      for c in axes(coeffs,2)
        coeff += coeffs[ldof,c]*proj[c,b]
      end
      aggdof_to_coeffs_data[p+b] = coeff
    end
  end
  nothing
end

function _fill_aggdof_to_coeffs_data_v1!(
  aggdof_to_coeffs_data,
  aggdof_to_dofs_ptrs,
  aggdof_to_fdof,
  fdof_to_acell,
  fdof_to_ldof,
  acell_to_coeffs,
  acell_to_proj,
  acell_to_is_owned,
  perm = eachindex(aggdof_to_fdof)
)
  cache2 = array_cache(acell_to_coeffs)
  cache3 = array_cache(acell_to_proj)
  for aggdof in perm
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    ! acell_to_is_owned[acell] && continue
    coeffs = getindex!(cache2,acell_to_coeffs,acell)
    proj = getindex!(cache3,acell_to_proj,acell)
    ldof = fdof_to_ldof[fdof]
    p = aggdof_to_dofs_ptrs[aggdof]-1
    row = @view coeffs[ldof,:]
    aggdof_to_coeffs_data[(p+1):(p+size(proj,2))] .= (row' * proj)'
  end
  nothing
end

function _fill_aggdof_to_coeffs_data_v2!(
  aggdof_to_coeffs_data,
  aggdof_to_dofs_ptrs,
  aggdof_to_fdof,
  fdof_to_acell,
  fdof_to_ldof,
  acell_to_coeffs,
  acell_to_proj,
  acell_to_is_owned,
  perm = eachindex(aggdof_to_fdof)
)
  cache2 = array_cache(acell_to_coeffs)
  cache3 = array_cache(acell_to_proj)
  for aggdof in perm
    fdof = aggdof_to_fdof[aggdof]
    acell = fdof_to_acell[fdof]
    ! acell_to_is_owned[acell] && continue
    ldof = fdof_to_ldof[fdof]
    p = aggdof_to_dofs_ptrs[aggdof]-1

    coeffs = getindex!(cache2,acell_to_coeffs,acell)
    proj = getindex!(cache3,acell_to_proj,acell)

    res = @view aggdof_to_coeffs_data[(p+1):(p+size(proj,2))]
    row = @view coeffs[ldof,:]
    #mul!(res, proj', row, 1.0, 0.0)
    mul!(res', row', proj, 1.0, 0.0)
  end
  nothing
end

function _fill_aggdof_to_coeffs_data_v3!(
  aggdof_to_coeffs_data,
  aggdof_to_dofs_ptrs,
  acells,
  acell_to_aggdofs,
  acell_to_coeffs,
  acell_to_proj,
  acell_to_is_owned
)
  cache2 = array_cache(acell_to_coeffs)
  cache3 = array_cache(acell_to_proj)
  for acell in acells
    ! acell_to_is_owned[acell] && continue
    coeffs = getindex!(cache2,acell_to_coeffs,acell)
    proj = getindex!(cache3,acell_to_proj,acell)
    agg_dofs = dataview(acell_to_aggdofs,acell)
    for (ldof,aggdof) in enumerate(agg_dofs)
      iszero(aggdof) && continue
      p = aggdof_to_dofs_ptrs[aggdof]-1
      row = @view coeffs[ldof,:]
      aggdof_to_coeffs_data[(p+1):(p+size(proj,2))] .= (row' * proj)'
    end
  end
  nothing
end

function _fill_aggdof_to_coeffs_data_v4!(
  aggdof_to_coeffs_data,
  aggdof_to_dofs_ptrs,
  acells,
  acell_to_aggdofs,
  acell_to_coeffs,
  acell_to_proj,
  acell_to_is_owned
)
  T = eltype(eltype(acell_to_coeffs))
  cache1 = CachedMatrix(T)
  cache2 = array_cache(acell_to_coeffs)
  cache3 = array_cache(acell_to_proj)
  for acell in acells
    ! acell_to_is_owned[acell] && continue
    coeffs = getindex!(cache2,acell_to_coeffs,acell)
    proj = getindex!(cache3,acell_to_proj,acell)
    setsize!(cache1, (size(proj,2), size(coeffs,2)))
    res = cache1.array
    mul!(res, proj', coeffs, 1.0, 0.0)
    agg_dofs = dataview(acell_to_aggdofs,acell)
    for (ldof,aggdof) in enumerate(agg_dofs)
      iszero(aggdof) && continue
      p = aggdof_to_dofs_ptrs[aggdof]-1
      aggdof_to_coeffs_data[(p+1):(p+size(proj,2))] .= (@view res[:,ldof])
    end
  end
  nothing
end

function _fill_aggdof_to_coeffs_data_v5!(
  aggdof_to_coeffs_data,
  aggdof_to_dofs_ptrs,
  acells,
  acell_to_ldofs,
  acell_to_coeffs,
  acell_to_proj,
  acell_to_is_owned
)
  T = eltype(eltype(acell_to_coeffs))
  cache1 = CachedMatrix(T)
  cache2 = array_cache(acell_to_coeffs)
  cache3 = array_cache(acell_to_proj)
  p = 1
  for acell in acells
    ! acell_to_is_owned[acell] && continue
    coeffs = getindex!(cache2,acell_to_coeffs,acell)
    proj = getindex!(cache3,acell_to_proj,acell)
    ldofs = dataview(acell_to_ldofs,acell)
    nl = length(ldofs)

    setsize!(cache1, (size(proj,2), nl))
    res = cache1.array
    mul!(res', (@view coeffs[ldofs,:]), proj, 1.0, 0.0)

    # WARNING: This relies on the fact that the aggdofs 
    # corresponding to a cell are contiguous in memory !!!
    range = aggdof_to_dofs_ptrs[p]:(aggdof_to_dofs_ptrs[p+nl]-1)
    @inbounds aggdof_to_coeffs_data[range] .= vec(res)

    p += nl
  end
  nothing
end

n = 20
geom = disk(1.4,x0=Point(1.0,1.7))
bgmodel = CartesianDiscreteModel((0,1,0,1),(n,n))
cutgeo = cut(bgmodel,geom)

trian = Triangulation(cutgeo,ACTIVE)

order = 3
reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(trian,reffe)

data = setup_aggdata_v0(cutgeo,V,V)

_fill_aggdof_to_coeffs_data_v0!(
  data.aggdof_to_coeffs.data, data.aggdof_to_dofs.ptrs, data.aggdof_to_fdof, data.fdof_to_acell, data.fdof_to_ldof, 
  data.acell_to_coeffs, data.acell_to_proj, data.acell_to_is_owned
)

cellsort(fdof) = (data.fdof_to_acell[fdof], fdof)
perm = sortperm(data.aggdof_to_fdof; by = cellsort)
_fill_aggdof_to_coeffs_data_v0!(
  data.aggdof_to_coeffs.data, data.aggdof_to_dofs.ptrs, data.aggdof_to_fdof, data.fdof_to_acell, data.fdof_to_ldof, 
  data.acell_to_coeffs, data.acell_to_proj, data.acell_to_is_owned, perm
)

# Optimizing cache hits by sorting aggdofs
function bm_v0(data,perm)
  _fill_aggdof_to_coeffs_data_v0!(
    data.aggdof_to_coeffs.data, data.aggdof_to_dofs.ptrs, data.aggdof_to_fdof, data.fdof_to_acell, data.fdof_to_ldof, 
    data.acell_to_coeffs, data.acell_to_proj, data.acell_to_is_owned, perm
  )
end
bm_v0(data, perm)
@benchmark bm_v0($data, $(eachindex(data.aggdof_to_fdof)))
@benchmark bm_v0($data, $perm)

# Optimizing vectorization by using views and BLAS
function bm_v1(data,perm)
  _fill_aggdof_to_coeffs_data_v1!(
    data.aggdof_to_coeffs.data, data.aggdof_to_dofs.ptrs, data.aggdof_to_fdof, data.fdof_to_acell, data.fdof_to_ldof, 
    data.acell_to_coeffs, data.acell_to_proj, data.acell_to_is_owned, perm
  )
end
bm_v1(data, perm)
@benchmark bm_v1($data, $perm)

function bm_v2(data,perm)
  _fill_aggdof_to_coeffs_data_v2!(
    data.aggdof_to_coeffs.data, data.aggdof_to_dofs.ptrs, data.aggdof_to_fdof, data.fdof_to_acell, data.fdof_to_ldof, 
    data.acell_to_coeffs, data.acell_to_proj, data.acell_to_is_owned, perm
  )
end
bm_v2(data, perm)
@benchmark bm_v2($data, $perm)

# Iterating over cells to avoid dictionary lookups when accessing the same array multiple times
fdof_to_aggdof = find_inverse_index_map(data.aggdof_to_fdof, length(data.fdof_to_acell))
aggdof_to_acell = data.fdof_to_acell[data.aggdof_to_fdof]
acells = sort!(unique(aggdof_to_acell))
is_aggdof = falses(length(data.fdof_to_acell)); is_aggdof[data.aggdof_to_fdof] .= true;
acell_to_dofs = get_cell_dof_ids(V)
acell_to_aggdofs = Table(collect(Int32,fdof_to_aggdof[acell_to_dofs.data]),acell_to_dofs.ptrs)

function bm_v3(data,acells,acell_to_aggdofs)
  _fill_aggdof_to_coeffs_data_v3!(
    data.aggdof_to_coeffs.data, data.aggdof_to_dofs.ptrs, acells, acell_to_aggdofs, 
    data.acell_to_coeffs, data.acell_to_proj, data.acell_to_is_owned
  )
end
bm_v3(data, acells, acell_to_aggdofs)
@benchmark bm_v3($data, $acells, $acell_to_aggdofs)

# Computing the whole matrix-matrix product. This will get less worth it 
# as the number of aggdofs per slave cell goes down.
function bm_v4(data,acells,acell_to_aggdofs)
  _fill_aggdof_to_coeffs_data_v4!(
    data.aggdof_to_coeffs.data, data.aggdof_to_dofs.ptrs, acells, acell_to_aggdofs, 
    data.acell_to_coeffs, data.acell_to_proj, data.acell_to_is_owned
  )
end
bm_v4(data, acells, acell_to_aggdofs)
@benchmark bm_v4($data, $acells, $acell_to_aggdofs)

# Lets actually reorder the aggdofs in memory to improve locality.

old_aggdof_to_acell = data.fdof_to_acell[data.aggdof_to_fdof]
new_aggdof_to_old_aggdof = sortperm(old_aggdof_to_acell)

data = setup_aggdata_v0(cutgeo,V,V)
data = (;
  aggdof_to_fdof = data.aggdof_to_fdof[new_aggdof_to_old_aggdof],
  aggdof_to_dofs = data.aggdof_to_dofs[new_aggdof_to_old_aggdof],
  aggdof_to_coeffs = data.aggdof_to_coeffs[new_aggdof_to_old_aggdof],
  fdof_to_acell = data.fdof_to_acell,
  fdof_to_ldof = data.fdof_to_ldof,
  acell_to_gcell = data.acell_to_gcell,
  acell_to_coeffs = data.acell_to_coeffs,
  acell_to_proj = data.acell_to_proj,
  acell_to_is_owned = data.acell_to_is_owned
)

fdof_to_aggdof = find_inverse_index_map(data.aggdof_to_fdof, length(data.fdof_to_acell))
aggdof_to_acell = data.fdof_to_acell[data.aggdof_to_fdof]
acells = sort!(unique(aggdof_to_acell))
is_aggdof = falses(length(data.fdof_to_acell)); is_aggdof[data.aggdof_to_fdof] .= true;
acell_to_dofs = get_cell_dof_ids(V)

acell_to_aggdofs = map(enumerate(acell_to_dofs)) do (acell,fdofs)
  aggdofs = map(fdofs) do fdof
    data.fdof_to_acell[fdof] == acell ? fdof_to_aggdof[fdof] : 0
  end
  collect(Int32,aggdofs)
end |> Table
@assert all(map(x -> issorted(findall(!iszero,x)), acell_to_aggdofs))

acell_to_ldofs = Table(map(x -> findall(!iszero,x), acell_to_aggdofs))

function bm_v5(data,acells,acell_to_ldofs)
  _fill_aggdof_to_coeffs_data_v5!(
    data.aggdof_to_coeffs.data, data.aggdof_to_dofs.ptrs, acells, acell_to_ldofs, 
    data.acell_to_coeffs, data.acell_to_proj, data.acell_to_is_owned
  )
end
bm_v5(data, acells, acell_to_ldofs)
@benchmark bm_v5($data, $acells, $acell_to_ldofs)

############################################################################################

dofs = get_fe_dof_basis(V)
basis = get_fe_basis(V)

Vg = V
dofs_g = get_fe_dof_basis(Vg)
basis_g = get_fe_basis(Vg)
phys_basis_g = change_domain(basis_g,PhysicalDomain())

cell_to_root = [7,6,6,7,6,6,7]

root_basis_g = GenericCellField(
  lazy_map(Reindex(get_data(phys_basis_g)),cell_to_root),trian,PhysicalDomain()
)

cell_to_coeffs = dofs(root_basis_g)
cell_to_proj = dofs_g(basis)

print_op_tree(phys_basis_g.cell_basis)

get_cell_map(get_grid(bgmodel))

