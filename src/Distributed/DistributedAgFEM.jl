
function AgFEMSpace(
  bgmodel::DistributedDiscreteModel,
  f::DistributedFESpace,
  bgcell_to_bgcellin::AbstractArray{<:AbstractVector},
  g::DistributedFESpace=f)

  bgmodel_gids = get_cell_gids(bgmodel)
  spaces = map(
    local_views(f),
    bgcell_to_bgcellin,
    local_views(g),
    local_views(bgmodel_gids)) do f,bgcell_to_bgcellin,g,gids
      AgFEMSpace(f,bgcell_to_bgcellin,g,local_to_global(gids))
  end
  trians = map(get_triangulation,local_views(f))
  trian = DistributedTriangulation(trians,bgmodel)
  trian = add_ghost_cells(trian)
  trian_gids = generate_cell_gids(trian)
  cell_to_cellin = _active_aggregates(bgcell_to_bgcellin)
  cell_to_ldofs = cell_ldof_to_mdof(spaces,cell_to_cellin)
  nldofs = map(num_free_dofs,spaces)
  gids = generate_gids(trian_gids,cell_to_ldofs,nldofs)
  vector_type = _find_vector_type(spaces,gids)
  DistributedSingleFieldFESpace(spaces,gids,trian,vector_type)
end

function aggregate(strategy,cutgeo::DistributedEmbeddedDiscretization,args...)
  aggregates,aggregate_owner = distributed_aggregate(strategy,cutgeo,args...)
  bgmodel = get_background_model(cutgeo)
  if has_remote_aggregation(bgmodel,aggregates)
    bgmodel = add_remote_aggregates(bgmodel,aggregates,aggregate_owner)
    cutgeo = change_bgmodel(cutgeo,bgmodel)
    aggregates = change_bgmodel(aggregates,get_cell_gids(bgmodel))
  end
  laggregates = _local_aggregates(aggregates,get_cell_gids(bgmodel))
  bgmodel,cutgeo,laggregates
end

function distributed_aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::DistributedEmbeddedDiscretization,
  in_or_out=IN)

  geo = get_geometry(cut)
  distributed_aggregate(strategy,cut,geo,in_or_out)
end

function distributed_aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out=IN)

  bgmodel = get_background_model(cut)
  facet_to_inoutcut = compute_bgfacet_to_inoutcut(bgmodel,geo)
  _distributed_aggregate_by_threshold(strategy.threshold,cut,geo,in_or_out,facet_to_inoutcut)
end


function _distributed_aggregate_by_threshold(threshold,cutgeo,geo,loc,facet_to_inoutcut)
  @assert loc in (IN,OUT)

  cutinorout = loc == IN ? (CUT_IN,IN) : (CUT_OUT,OUT)
  trian = Triangulation(cutgeo,cutinorout,geo)
  model = get_background_model(cutgeo)
  bgtrian = get_triangulation(model)
  cell_to_cut_meas = map(_get_cell_measure,local_views(trian),local_views(bgtrian))
  cell_to_meas = map(get_cell_measure,local_views(bgtrian))
  cell_to_unit_cut_meas = map(cell_to_cut_meas,cell_to_meas) do c_to_cm,c_to_m
    lazy_map(/,c_to_cm,c_to_m)
  end

  cell_to_inoutcut = compute_bgcell_to_inoutcut(cutgeo,geo)

  cell_to_coords = map(get_cell_coordinates,local_views(model))
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  cell_to_faces = map(t->get_faces(t,D,D-1),local_views(topo))
  face_to_cells = map(t->get_faces(t,D-1,D),local_views(topo))
  gids = get_cell_gids(model)

  _distributed_aggregate_by_threshold_barrier(
    threshold,cell_to_unit_cut_meas,facet_to_inoutcut,cell_to_inoutcut,
    loc,cell_to_coords,cell_to_faces,face_to_cells,gids)
end


function _distributed_aggregate_by_threshold_barrier(
  threshold,cell_to_unit_cut_meas,facet_to_inoutcut,cell_to_inoutcut,
  loc,cell_to_coords,cell_to_faces,face_to_cells,gids)


  ocell_to_touched = map(cell_to_unit_cut_meas) do c_to_m
    map(≥,c_to_m,Fill(threshold,length(c_to_m)))
  end
  cell_to_touched = _add_ghost_values(ocell_to_touched,gids)

  cell_to_root_centroid = map(cell_to_coords) do cell_to_coords
    map(i->sum(i)/length(i),cell_to_coords)
  end
  PVector(cell_to_root_centroid,partition(gids)) |> consistent! |> wait

  n_cells = map(length,cell_to_touched)
  touched_cells = map(findall,cell_to_touched)

  cell_to_cellin = map(n->zeros(Int32,n),n_cells)
  map(cell_to_cellin,touched_cells,local_to_global(gids)) do c_to_ci,cells,l_to_g
    gcells = lazy_map(Reindex(l_to_g),cells)
    c_to_ci[cells] = gcells
  end

  cell_to_neig = map(n->zeros(Int32,n),n_cells)

  cell_to_root_part = map(collect,local_to_owner(gids))

  c1 = map(array_cache,cell_to_faces)
  c2 = map(array_cache,face_to_cells)

  max_iters = 20
  for iter in 1:max_iters
    all_aggregated = _aggregate_one_step!(c1,c2,gids,
      cell_to_inoutcut,
      cell_to_touched,
      cell_to_neig,
      cell_to_cellin,
      cell_to_root_centroid,
      cell_to_root_part,
      cell_to_faces,
      face_to_cells,
      facet_to_inoutcut,
      loc)

    PVector(cell_to_touched,partition(gids)) |> consistent! |> wait
    PVector(cell_to_neig,partition(gids)) |> consistent! |> wait
    PVector(cell_to_cellin,partition(gids)) |> consistent! |> wait
    PVector(cell_to_root_centroid,partition(gids)) |> consistent! |> wait
    PVector(cell_to_root_part,partition(gids)) |> consistent! |> wait

    reduction!(&,all_aggregated,all_aggregated,destination=:all)

    if PartitionedArrays.getany(all_aggregated)
      break
    end
  end

  cell_to_cellin, cell_to_root_part, cell_to_neig
end

function _aggregate_one_step!(c1,c2,gids::PRange,
  cell_to_inoutcut,
  cell_to_touched,
  cell_to_neig,
  cell_to_cellin,
  cell_to_root_centroid,
  cell_to_root_part,
  cell_to_faces,
  face_to_cells,
  facet_to_inoutcut,
  loc)

  map(c1,c2,own_to_local(gids),
    cell_to_inoutcut,
    cell_to_touched,
    cell_to_neig,
    cell_to_cellin,
    cell_to_root_centroid,
    cell_to_root_part,
    local_to_global(gids),
    cell_to_faces,
    face_to_cells,
    facet_to_inoutcut) do c1,c2,own_cells,
        cell_to_inoutcut,
        cell_to_touched,
        cell_to_neig,
        cell_to_cellin,
        cell_to_root_centroid,
        cell_to_root_part,
        cell_to_gcell,
        cell_to_faces,
        face_to_cells,
        facet_to_inoutcut

    _aggregate_one_step!(
      c1,c2,own_cells,
      cell_to_inoutcut,
      cell_to_touched,
      cell_to_neig,
      cell_to_cellin,
      cell_to_root_centroid,
      cell_to_root_part,
      cell_to_gcell,
      cell_to_faces,
      face_to_cells,
      facet_to_inoutcut,
      loc)
  end
end

function _aggregate_one_step!(
  c1,c2,own_cells,
  cell_to_inoutcut,
  cell_to_touched,
  cell_to_neig,
  cell_to_cellin,
  cell_to_root_centroid,
  cell_to_root_part,
  cell_to_gcell,
  cell_to_faces,
  face_to_cells,
  facet_to_inoutcut,
  loc)

  all_aggregated = true
  for cell in own_cells
    if ! cell_to_touched[cell] && cell_to_inoutcut[cell] == CUT
      neigh_cell = _find_best_neighbor_from_centroid_distance(
        c1,c2,cell,
        cell_to_faces,
        face_to_cells,
        cell_to_touched,
        cell_to_root_centroid,
        facet_to_inoutcut,
        loc)
      if neigh_cell > 0
        cellin = cell_to_cellin[neigh_cell]
        centroid = cell_to_root_centroid[neigh_cell]
        part = cell_to_root_part[neigh_cell]
        neigh_gcell = cell_to_gcell[neigh_cell]

        cell_to_neig[cell] = neigh_gcell
        cell_to_cellin[cell] = cellin
        cell_to_root_centroid[cell] = centroid
        cell_to_root_part[cell] = part
      else
        all_aggregated = false
      end
    end
  end
  _touch_aggregated_cells!(cell_to_touched,cell_to_cellin)
  all_aggregated
end

function _find_best_neighbor_from_centroid_distance(
  c1,c2,cell,
  cell_to_faces,
  face_to_cells,
  cell_to_touched,
  cell_to_root_centroid,
  facet_to_inoutcut,
  loc)

  faces = getindex!(c1,cell_to_faces,cell)
  dmin = Inf
  T = eltype(eltype(face_to_cells))
  best_neigh_cell = zero(T)
  for face in faces
    inoutcut = facet_to_inoutcut[face]
    if  inoutcut != CUT && inoutcut != loc
      continue
    end
    neigh_cells = getindex!(c2,face_to_cells,face)
    for neigh_cell in neigh_cells
      if neigh_cell != cell && cell_to_touched[neigh_cell]
        p = cell_to_root_centroid[neigh_cell]
        q = cell_to_root_centroid[cell]
        d = norm(p-q)
        if (1.0+1.0e-9)*d < dmin
          dmin = d
          best_neigh_cell = neigh_cell
        end
      end
    end
  end
  best_neigh_cell
end

function _add_ghost_values(own_v,gids::PRange)
  lens = map(length,local_views(gids))
  eltypes = map(eltype,own_v)
  local_v = map(zeros,eltypes,lens)
  map(local_v,own_v,own_to_local(gids)) do l,o,o_to_l
    l[o_to_l] = o
  end
  PVector(local_v,partition(gids)) |> consistent! |> wait
  local_v
end

function _get_cell_measure(trian1::Triangulation,trian2::Triangulation)
  if num_cells(trian1) == 0
    Fill(0.0,num_cells(trian2))
  else
    get_cell_measure(trian1,trian2)
  end
end

function merge_nodes(model::DistributedDiscreteModel)
  node_gids = get_face_gids(model,0)
  cell_gids = get_cell_gids(model)
  models = map(local_views(model),local_views(node_gids)) do model,node_ids
    merge_nodes(model,node_ids)
  end
  DistributedDiscreteModel(models,cell_gids)
end

function merge_nodes(model::DiscreteModel,ids)
  l_to_g = local_to_global(ids)
  n_global = length(global_to_local(ids))
  n_local = length(l_to_g)
  g_to_l = VectorFromDict(reverse(l_to_g),reverse(1:length(l_to_g)),n_global)
  l_to_lparent = g_to_l[l_to_g]
  lnew_to_l = findall(map(==,l_to_lparent,1:n_local))
  l_to_lnew = zeros(Int,n_local)
  l_to_lnew[lnew_to_l] = 1:length(lnew_to_l)
  g_to_lnew = map(l->iszero(l) ? l : l_to_lnew[l], g_to_l)
  l_to_lnew = g_to_lnew[l_to_g]

  grid = get_grid(model)
  coords = get_node_coordinates(grid)
  conn = get_cell_node_ids(grid)
  reffes = get_reffes(grid)
  ctypes = get_cell_type(grid)

  coords = coords[lnew_to_l]
  conn = map(Broadcasting(Reindex(l_to_lnew)),conn) |> Table
  grid = UnstructuredGrid(coords,conn,reffes,ctypes)
  UnstructuredDiscreteModel(grid)
end

function add_remote_cells(model::DistributedDiscreteModel,remote_cells,remote_parts)
  # Send remote gids to owners
  snd_ids = remote_parts
  snd_remotes = remote_cells
  graph = ExchangeGraph(snd_ids)
  rcv_remotes = allocate_exchange(snd_remotes,graph)
  exchange!(rcv_remotes,snd_remotes,graph) |> wait

  # Send remote coordinates
  gids = get_cell_gids(model)
  snd_gids = rcv_remotes
  snd_lids = map(global_to_local(gids),snd_gids) do g_to_l,gids
    map(Reindex(g_to_l),gids)
  end
  snd_coords = map(local_views(model),snd_lids) do m,lids
    T = eltype(eltype(get_cell_coordinates(m)))
    coords = map(lids) do lids
      coords = map(Reindex(get_cell_coordinates(m)),lids)
      reduce(append!,coords,init=T[])
    end
    Vector{Vector{T}}(coords)
  end
  rgraph = reverse(graph)
  rcv_coords = allocate_exchange(snd_coords,rgraph)
  exchange!(rcv_coords,snd_coords,rgraph) |> wait

  # Build remote grids
  ncells = map(remote_cells) do cells
    sum(length,cells,init=0)
  end
  reffes = map(get_reffes,local_views(model))
  reffe = map(only,reffes)
  ctypes = map(n->ones(Int,n),ncells)
  coords = map(PartitionedArrays.getdata,rcv_coords)
  conn = map(ncells,reffe) do ncells,reffe
    n = num_nodes(reffe)
    data = 1:n*ncells
    ptrs = 1:n:n*ncells+1
    Table(data,ptrs)
  end
  rgrids = map(UnstructuredGrid,coords,conn,reffes,ctypes)

  # Build appended model
  lgrids = map(get_grid,local_views(model))
  _grids = map(lazy_append,lgrids,rgrids)
  _models = map(UnstructuredDiscreteModel,_grids)
  agids = add_remote_ids(gids,remote_cells,remote_parts)
  amodel = DistributedDiscreteModel(_models,agids) |> merge_nodes
<<<<<<< HEAD
  
=======

>>>>>>> d579786312ab8726581d7aa664f3a61c4d1c2b11
  D = num_cell_dims(model)
  grids = map(get_grid,local_views(amodel))
  topos = map(get_grid_topology,local_views(amodel))
  d_to_dface_to_entity = map(topos) do topo
    [ fill(Int32(UNSET),num_faces(topo,d)) for d in 0:D ]
  end
  oldtopos = map(get_grid_topology,local_views(model))
  oldlabels = map(get_face_labeling,local_views(model))
  for d in 0:D
    _fill_labels!(d_to_dface_to_entity,
                  oldlabels,
                  oldtopos,
                  topos,
                  ncells,
                  d,D,
                  snd_lids,
                  rgraph)
  end
  labels = map(d_to_dface_to_entity,oldlabels) do d_to_dface_to_entity,ol
    FaceLabeling(d_to_dface_to_entity,ol.tag_to_entities,ol.tag_to_name)
  end
  models = map(UnstructuredDiscreteModel,grids,topos,labels)
  DistributedDiscreteModel(models,agids)
end

function _fill_labels!(d_to_dface_to_entity,llabels,ltopos,ntopos,nremotes,
                       d::Int,D::Int,snd_lids,rgraph)
  snd_labels = map(ltopos,llabels,snd_lids) do topo,labels,lids
    dface_to_entity = get_face_entity(labels,d)
    T = eltype(eltype(dface_to_entity))
    labels = map(lids) do lids
      faces = map(Reindex(get_faces(topo,D,d)),lids)
      labels = map(Reindex(dface_to_entity),faces)
      reduce(append!,labels,init=T[])
    end
    Vector{Vector{T}}(labels)
  end
  rcv_labels = allocate_exchange(snd_labels,rgraph)
  exchange!(rcv_labels,snd_labels,rgraph) |> wait
  rlabels = map(PartitionedArrays.getdata,rcv_labels)
  map(d_to_dface_to_entity,llabels) do d_to_dface_to_entity,l
    l_face_entity = get_face_entity(l,d)
    d_to_dface_to_entity[d+1][1:length(l_face_entity)] = l_face_entity
  end
  map(d_to_dface_to_entity,ntopos,nremotes,rlabels) do d_to_dface_to_entity,nt,nr,r
    nr > 0 && begin
      rf = reduce(vcat,get_faces(nt,D,d)[end-nr+1:end])
      d_to_dface_to_entity[d+1][rf] = r
    end
  end
end

function add_remote_aggregates(model::DistributedDiscreteModel,aggregates,aggregate_owner)
  gids = get_cell_gids(model)
  remote_cells,remote_parts = _extract_remote_cells(gids,aggregates,aggregate_owner)
  remote_cells,remote_parts = _group_remote_ids(remote_cells,remote_parts)
  add_remote_cells(model,remote_cells,remote_parts)
end

function _extract_remote_cells(gids::PRange,aggregates,aggregate_owner)
  remote_aggids = map(aggregates,global_to_local(gids)) do agg,g_to_l
    ids = findall(agg) do i
      !iszero(i) && iszero(g_to_l[i])
    end
    unique(Reindex(agg),ids)
  end

  remote_cells = map(aggregates,remote_aggids) do agg,ids
    map(Reindex(agg),ids)
  end

  remote_parts = map(aggregate_owner,remote_aggids) do agg,ids
    map(Reindex(agg),ids)
  end

  remote_cells,remote_parts
end

function _group_remote_ids(remote_ids,remote_parts)
  new_parts = map(sort∘unique,remote_parts)
  new_ids = map(remote_ids,remote_parts,new_parts) do ids,parts,uparts
    grouped_ids = map(i->Int[],1:length(uparts))
    for (id,p) in zip(ids,parts)
      j = findfirst(==(p),uparts)
      union!(grouped_ids[j],id)
    end
    map!(sort!,grouped_ids,grouped_ids)
  end
  new_ids,new_parts
end

function _ungroup_remote_ids(remote_ids,remote_parts)
  new_ids = map(remote_ids) do ids
    reduce(append!,ids,init=eltype(eltype(ids))[])
  end
  new_parts = map(remote_ids,remote_parts) do ids,parts
    n = map(length,ids)
    parts_v = map((p,n)->Fill(p,n),parts,n)
    reduce(append!,parts_v,init=eltype(parts)[])
  end
  new_ids,new_parts
end

function add_remote_ids(gids::PRange,remote_gids,remote_parts)
  new_gids,new_parts = _ungroup_remote_ids(remote_gids,remote_parts)
  lid_to_gid = map(vcat,local_to_global(gids),new_gids)
  lid_to_part = map(vcat,local_to_owner(gids),new_parts)
  p = map(lid_to_gid,lid_to_part,partition(gids)) do l_to_g,l_to_p,p
    l_to_g = collect(Int,l_to_g)
    l_to_p = collect(Int32,l_to_p)
    LocalIndices(length(gids),part_id(p),l_to_g,l_to_p)
  end
  PRange(p)
end

function has_remote_aggregation(model::DistributedDiscreteModel,aggregates)
  gids = get_cell_gids(model)
  has_remote_aggregation(aggregates,gids)
end

function has_remote_aggregation(aggregates,gids::PRange)
  remote_aggregation = map(aggregates,global_to_local(gids)) do agg,g_to_l
    lazy_map(agg) do a
    iszero(a) || !iszero(g_to_l[a])
    end |> all |> !
  end
  reduction(|,remote_aggregation,destination=:all) |> PartitionedArrays.getany
end


function _active_aggregates(bgcell_to_bgcellin::AbstractVector{<:AbstractVector})
  map(_active_aggregates,bgcell_to_bgcellin)
end

function _active_aggregates(bgcell_to_bgcellin)
  acell_to_bgcell = findall(!iszero,bgcell_to_bgcellin)
  bgcell_to_acell = zeros(Int,length(bgcell_to_bgcellin))
  bgcell_to_acell[acell_to_bgcell] = 1:length(acell_to_bgcell)
  acell_to_bgcellin = bgcell_to_bgcellin[ acell_to_bgcell ]
  bgcell_to_acell[ acell_to_bgcellin ]
end

function cell_ldof_to_mdof(
  spaces::AbstractArray{<:FESpace},
  cell_to_cellin::AbstractArray{<:AbstractVector})

  map(cell_ldof_to_mdof,spaces,cell_to_cellin)
end

function cell_ldof_to_mdof(
  space::FESpaceWithLinearConstraints,
  cell_to_cellin::AbstractVector)

  DOF_to_mDOFs = space.DOF_to_mDOFs
  cell_ldof_to_dof = space.cell_to_ldof_to_dof
  n_fdofs = space.n_fdofs
  n_fmdofs = space.n_fmdofs
  cell_ldof_to_mdof = map(cell_ldof_to_dof) do ldof_to_dof
    map(ldof_to_dof) do dof
      DOF = _dof_to_DOF(dof,n_fdofs)
      mDOFs = DOF_to_mDOFs[DOF]
      length(mDOFs) == 1 ? _DOF_to_dof(mDOFs[1],n_fmdofs) : zero(eltype(mDOFs))
    end
  end
  for (cell,ldof_to_mdof) in enumerate(cell_ldof_to_mdof)
    if cell_to_cellin[cell] != cell
     empty!(ldof_to_mdof)
    end
  end
  cell_ldof_to_mdof
end

function _local_aggregates(cell_to_gcellin,gids::PRange)
  map(_local_aggregates,cell_to_gcellin,global_to_local(gids))
end

function _local_aggregates(cell_to_gcellin,gcell_to_cell)
  T = eltype(cell_to_gcellin)
  map(cell_to_gcellin) do gcin
    iszero(gcin) ? gcin : T(gcell_to_cell[ gcin ])
  end
end

function change_bgmodel(cell_to_gcellin,gids::PRange)
  map(change_bgmodel,cell_to_gcellin,local_to_global(gids))
end

function change_bgmodel(cell_to_gcellin,ncell_to_gcell)
  ncells = length(cell_to_gcellin)
  ncell_to_gcellin = zeros(Int,length(ncell_to_gcell))
  for (ncell,gcell) in enumerate(ncell_to_gcell)
    ncell_to_gcellin[ncell] = ncell > ncells ? gcell : cell_to_gcellin[ncell]
  end
  ncell_to_gcellin
end

# Specialised methods for Algoim quadratures

function aggregate(bgmodel::DistributedDiscreteModel,
                   cell_to_is_active,
                   cell_to_is_cut,
                   in_or_out)
  aggregates,aggregate_owner = distributed_aggregate(bgmodel,
                                                     cell_to_is_active,
                                                     cell_to_is_cut,
                                                     in_or_out)
  if has_remote_aggregation(bgmodel,aggregates)
    bgmodel = add_remote_aggregates(bgmodel,aggregates,aggregate_owner)
    aggregates = change_bgmodel(aggregates,get_cell_gids(bgmodel))
  end
  laggregates = _local_aggregates(aggregates,get_cell_gids(bgmodel))
  bgmodel,laggregates
end

function aggregate(bgtrian::DistributedTriangulation,
                   cell_to_is_active,
                   cell_to_is_cut,
                   in_or_out)
  msg = """
  This is not implemented, because the aggregation might need
  to extend the background model with the remote aggregates.
  """
  @notimplemented  msg
end

function distributed_aggregate(bgmodel::DistributedDiscreteModel,
                               cell_to_is_active,
                               cell_to_is_cut,
                               in_or_out)
  n_cells = length(cell_to_is_active)
  @assert n_cells == length(cell_to_is_cut)

  cell_to_unit_cut_meas = map(cell_to_is_active,cell_to_is_cut) do visa, visc
    lazy_map(visa,visc) do isa, isc
      !isa ? 0.0 : (isc ? 0.0 : 1.0)
    end
  end

  ocell_to_inoutcut = map(cell_to_is_active,cell_to_is_cut) do visa, visc
    lazy_map(visa,visc) do isa, isc
      !isa ? OUT : (isc ? CUT : IN)
    end
  end

  cell_to_coords = map(get_cell_coordinates,local_views(bgmodel))
  topo = get_grid_topology(bgmodel)
  D = num_cell_dims(bgmodel)
  cell_to_faces = map(t->get_faces(t,D,D-1),local_views(topo))
  face_to_cells = map(t->get_faces(t,D-1,D),local_views(topo))
  gids = get_cell_gids(bgmodel)
  # A hack follows to avoid constructing the actual facet_to_inoutcut array
  facet_to_inoutcut = map(m->fill(in_or_out,num_faces(m,D-1)),local_views(bgmodel))

  threshold = 1.0
  cell_to_inoutcut = _add_ghost_values(ocell_to_inoutcut,gids)
  _distributed_aggregate_by_threshold_barrier(
    threshold,cell_to_unit_cut_meas,facet_to_inoutcut,cell_to_inoutcut,
    in_or_out,cell_to_coords,cell_to_faces,face_to_cells,gids)
end