
function AgFEMSpace(
  bgmodel::DistributedDiscreteModel,
  f::DistributedFESpace,
  bgcell_to_bgcellin::AbstractArray{<:AbstractVector},
  g::DistributedFESpace=f
)
  bgmodel_gids = get_cell_gids(bgmodel)
  spaces = map(
    local_views(f),
    bgcell_to_bgcellin,
    local_views(g),
    local_views(bgmodel_gids)) do f,bgcell_to_bgcellin,g,gids
      AgFEMSpace(f,bgcell_to_bgcellin,g,local_to_global(gids))
  end
  trian = add_ghost_cells(get_triangulation(f))
  trian_gids = generate_cell_gids(trian)
  cell_to_cellin = _active_aggregates(bgcell_to_bgcellin)
  cell_to_ldofs = cell_ldof_to_mdof(spaces,cell_to_cellin)
  nldofs = map(num_free_dofs,spaces)
  gids = generate_gids(trian_gids,cell_to_ldofs,nldofs)
  vector_type = _find_vector_type(spaces,gids)
  DistributedSingleFieldFESpace(spaces,gids,trian,vector_type)
end

function aggregate(strategy,cutgeo::DistributedEmbeddedDiscretization,args...)
  aggregates, aggregate_owner = distributed_aggregate(strategy,cutgeo,args...)
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
  in_or_out=IN
)
  geo = get_geometry(cut)
  distributed_aggregate(strategy,cut,geo,in_or_out)
end

function distributed_aggregate(
  strategy::AggregateCutCellsByThreshold,
  cut::DistributedEmbeddedDiscretization,
  geo::CSG.Geometry,
  in_or_out=IN
)
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
    loc,cell_to_coords,cell_to_faces,face_to_cells,gids
  )
end

function _distributed_aggregate_by_threshold_barrier(
  threshold,cell_to_unit_cut_meas,facet_to_inoutcut,cell_to_inoutcut,
  loc,cell_to_coords,cell_to_faces,face_to_cells,gids
)
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
      loc
    )

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
  loc
)
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

# change_bgmodel

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

function change_bgmodel(
  cutgeo::DistributedEmbeddedDiscretization,
  model::DistributedDiscreteModel
)
  cuts = map(change_bgmodel,local_views(cutgeo),local_views(model))
  DistributedEmbeddedDiscretization(cuts,model)
end

function change_bgmodel(
  cutgeo::DistributedEmbeddedDiscretization,
  model::DistributedDiscreteModel,
  cell_to_new_cell
)
  cuts = map(change_bgmodel,local_views(cutgeo),local_views(model),cell_to_new_cell)
  DistributedEmbeddedDiscretization(cuts,model)
end

function change_bgmodel(
  cut::EmbeddedDiscretization,
  newmodel::DiscreteModel,
  cell_to_newcell=1:num_cells(get_background_model(cut))
)
  ls_to_bgc_to_ioc = map(cut.ls_to_bgcell_to_inoutcut) do bgc_to_ioc
    new_bgc_to_ioc = Vector{Int8}(undef,num_cells(newmodel))
    new_bgc_to_ioc[cell_to_newcell] = bgc_to_ioc
    new_bgc_to_ioc
  end
  subcells = change_bgmodel(cut.subcells,cell_to_newcell)
  subfacets = change_bgmodel(cut.subfacets,cell_to_newcell)
  EmbeddedDiscretization(
    newmodel,
    ls_to_bgc_to_ioc,
    subcells,
    cut.ls_to_subcell_to_inout,
    subfacets,
    cut.ls_to_subfacet_to_inout,
    cut.oid_to_ls,
    cut.geo
  )
end

function change_bgmodel(
  cut::EmbeddedFacetDiscretization,
  newmodel::DiscreteModel,
  facet_to_newfacet=1:num_facets(get_background_model(cut))
)
  nfacets = num_facets(newmodel)
  ls_to_bgf_to_ioc = map(cut.ls_to_facet_to_inoutcut) do bgf_to_ioc
    new_bgf_to_ioc = Vector{Int8}(undef,nfacets)
    new_bgf_to_ioc[facet_to_newfacet] = bgf_to_ioc
    new_bgf_to_ioc
  end
  subfacets = change_bgmodel(cut.subfacets,facet_to_newfacet)
  EmbeddedFacetDiscretization(
    newmodel,
    ls_to_bgf_to_ioc,
    subfacets,
    cut.ls_to_subfacet_to_inout,
    cut.oid_to_ls,
    cut.geo
  )
end

function change_bgmodel(cells::SubCellData,cell_to_newcell)
  cell_to_bgcell = lazy_map(Reindex(cell_to_newcell),cells.cell_to_bgcell)
  SubCellData(
    cells.cell_to_points,
    collect(Int32,cell_to_bgcell),
    cells.point_to_coords,
    cells.point_to_rcoords
  )
end

function change_bgmodel(facets::SubFacetData,cell_to_newcell)
  facet_to_bgcell = lazy_map(Reindex(cell_to_newcell),facets.facet_to_bgcell)
  SubFacetData(
    facets.facet_to_points,
    facets.facet_to_normal,
    collect(Int32,facet_to_bgcell),
    facets.point_to_coords,
    facets.point_to_rcoords
  )
end

function AgFEMSpace(
  f::DistributedSingleFieldFESpace,
  bgcell_to_bgroot::AbstractVector,
  g::DistributedSingleFieldFESpace = f
)
  trian = get_triangulation(g)
  bgcell_to_gcell = map(local_views(trian)) do trian
    glue = get_glue(trian,Val(num_cell_dims(trian)))
    glue.mface_to_tface
  end
  AgFEMSpace(
    f,bgcell_to_bgroot,get_fe_basis(g),get_fe_dof_basis(g),bgcell_to_gcell
  )
end

function AgFEMSpace(
  f::DistributedSingleFieldFESpace,
  bgcell_to_bgroot::AbstractVector,
  shfns_g::DistributedCellField,
  dofs_g::DistributedCellDof,
  bgcell_to_gcell
)
  trian = get_triangulation(f)
  spaces = local_views(f) # Not even necessary to have a DistributedFESpace
  cell_gids = GridapDistributed.generate_cell_gids(trian)

  # The following is the same as in serial, it's just some reindexing onto the triangulation
  nlfdofs, cell_to_root, cell_to_fdofs, cell_to_coeffs, cell_to_proj, cell_to_gcell = map(
    spaces, bgcell_to_bgroot, local_views(shfns_g), local_views(dofs_g), bgcell_to_gcell
  ) do space, bgcell_to_bgroot, shfns_g, dofs_g, bgcell_to_gcell
    trian = get_triangulation(space)

    glue = get_glue(trian,Val(num_cell_dims(trian)))
    cell_to_bgcell = glue.tface_to_mface
    bgcell_to_cell = glue.mface_to_tface
    cell_to_bgroot = view(bgcell_to_bgroot,cell_to_bgcell)
    cell_to_root = collect(Int32,lazy_map(Reindex(bgcell_to_cell),cell_to_bgroot))
    cell_to_gcell = collect(Int32,lazy_map(Reindex(bgcell_to_gcell),cell_to_bgcell))

    cell_phys_shfns_g = get_array(change_domain(shfns_g,PhysicalDomain()))
    cell_phys_root_shfns_g = lazy_map(Reindex(cell_phys_shfns_g),cell_to_root)
    root_shfns_g = GenericCellField(cell_phys_root_shfns_g,trian,PhysicalDomain())

    # Compute data needed to compute the constraints
    dofs_f = get_fe_dof_basis(space)
    shfns_f = get_fe_basis(space)
    cell_to_coeffs = dofs_f(root_shfns_g)
    cell_to_proj = dofs_g(shfns_f)
    cell_to_fdofs = get_cell_dof_ids(space)
    nlfdofs = num_free_dofs(space)

    return nlfdofs, cell_to_root, cell_to_fdofs, cell_to_coeffs, cell_to_proj, cell_to_gcell
  end |> tuple_of_arrays

  fdof_is_agg, fdof_to_cell, fdof_to_ldof = map(
    nlfdofs, cell_to_root, cell_to_fdofs, cell_to_gcell
  ) do nlfdofs, cell_to_root, cell_to_fdofs, cell_to_gcell
    fdof_is_agg, fdof_to_cell, fdof_to_ldof = AgFEM._allocate_fdof_to_data(nlfdofs)
    AgFEM._fill_fdof_to_data!(fdof_is_agg,fdof_to_cell,fdof_to_ldof,cell_to_root,cell_to_fdofs,cell_to_gcell)
    return fdof_is_agg, fdof_to_cell, fdof_to_ldof
  end |> tuple_of_arrays

  mdof_gids, aggdof_gids, mdof_to_fdof, aggdof_to_fdof, fdof_to_posneg = generate_aggregated_gids(
    cell_gids, cell_to_fdofs, fdof_is_agg
  )

  # Create aggdof to nldofs mapping (i.e count how many masters every agg dof has)
  ptrs, aggdof_to_nldofs = map(
    aggdof_to_fdof,fdof_to_cell,cell_to_root,cell_to_fdofs,partition(aggdof_gids)
  ) do aggdof_to_fdof,fdof_to_cell,cell_to_root,cell_to_fdofs,aggdof_ids
    ptrs = get_aggdof_ptrs(
      aggdof_to_fdof,fdof_to_cell,cell_to_root,cell_to_fdofs,own_to_local(aggdof_ids)
    )
    aggdof_to_nldofs = view(ptrs,2:length(ptrs))
    return ptrs, aggdof_to_nldofs
  end |> tuple_of_arrays
  consistent!(PVector(aggdof_to_nldofs, partition(aggdof_gids))) |> wait

  # Fill the info for the owned aggdofs
  aggdof_to_dofs, aggdof_to_coeffs = map(
    ptrs,aggdof_to_fdof,fdof_to_cell,fdof_to_ldof,cell_to_root,cell_to_fdofs,cell_to_coeffs,cell_to_proj,partition(aggdof_gids)
  ) do ptrs,aggdof_to_fdof,fdof_to_cell,fdof_to_ldof,cell_to_root,cell_to_fdofs,cell_to_coeffs,cell_to_proj,aggdof_ids
    aggdof_to_dofs_data, aggdof_to_coeffs_data = AgFEM._allocate_aggdof_to_data(ptrs,cell_to_coeffs)
    aggdof_to_dofs!(
      aggdof_to_dofs_data,ptrs,aggdof_to_fdof,fdof_to_cell,cell_to_root,cell_to_fdofs,own_to_local(aggdof_ids)
    )
    aggdof_to_coeffs!(
      aggdof_to_coeffs_data,ptrs,aggdof_to_fdof,fdof_to_cell,fdof_to_ldof,cell_to_coeffs,cell_to_proj,own_to_local(aggdof_ids)
    )
    return JaggedArray(aggdof_to_dofs_data,ptrs), JaggedArray(aggdof_to_coeffs_data,ptrs)
  end |> tuple_of_arrays

  # Make consistent:
  #  - coeffs are straightforward
  #  - dofs have to be converted to mdof gids, communicated, then converted back to local dofs
  map(to_global_dofs!, aggdof_to_dofs, partition(mdof_gids), fdof_to_posneg)
  t1 = consistent!(PVector(aggdof_to_dofs, partition(aggdof_gids)))
  t2 = consistent!(PVector(aggdof_to_coeffs, partition(aggdof_gids)))
  wait(t1)
  expanded_dof_gids = PRange(
    map(to_local_dofs!,aggdof_to_dofs, partition(aggdof_gids), partition(mdof_gids), mdof_to_fdof, nlfdofs)
  )
  wait(t2)

  aggdof_to_dofs, aggdof_to_coeffs = map(aggdof_to_dofs, aggdof_to_coeffs) do dofs, coeffs
    Table(dofs.data,dofs.ptrs), Table(coeffs.data,coeffs.ptrs)
  end |> tuple_of_arrays
  # agg_spaces = map(FESpaceWithLinearConstraints, aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs, spaces)

  return aggdof_to_fdof, aggdof_to_dofs, aggdof_to_coeffs, aggdof_gids, expanded_dof_gids
end

# This is easy, all nonzero entries are local
function to_global_dofs!(aggdof_to_dofs, mdof_ids, dof_to_mdof)
  mdof_lid_to_gid = local_to_global(mdof_ids)
  data = aggdof_to_dofs.data
  for k in eachindex(data)
    iszero(data[k]) && continue
    mdof = dof_to_mdof[data[k]]
    data[k] = mdof_lid_to_gid[mdof]
  end
end

# This one is tricky: some nonzero entries will be non-local (i.e roots on other processors).
# We have to add these to the pre-existing dof numbering.
function to_local_dofs!(aggdof_to_dofs, aggdof_ids, mdof_ids, mdof_to_dof, n_dofs)
  rank = part_id(mdof_ids)
  n_mdofs = length(mdof_to_dof)
  new_gids = Dict{Int,Tuple{Int32,Int32}}()
  dof_gid_to_lid = global_to_local(mdof_ids)
  ptrs = aggdof_to_dofs.ptrs
  data = aggdof_to_dofs.data
  for (aggdof,owner) in enumerate(local_to_owner(aggdof_ids))
    for k in ptrs[aggdof]:ptrs[aggdof+1]-1
      @assert !iszero(data[k]) # All entries should be nonzero now
      gid = data[k]
      mdof = dof_gid_to_lid[gid]
      if !iszero(mdof) # Local dof
        dof = mdof_to_dof[mdof]
      else # Remote dof
        mdof, mdof_owner = get!(new_gids,gid,(n_mdofs+1,owner))
        dof = n_dofs + mdof - length(mdof_to_dof)
        @assert isequal(mdof_owner,owner) && !isequal(owner, rank)
        n_mdofs += isequal(mdof,n_mdofs+1) # Only increment if new
      end
      data[k] = dof
    end
  end
  
  # Create expanded master dof numbering
  lid_to_gid = Vector{Int}(undef,n_mdofs)
  lid_to_owner = Vector{Int32}(undef,n_mdofs)
  lid_to_gid[1:length(mdof_ids)] .= local_to_global(mdof_ids)
  lid_to_owner[1:length(mdof_ids)] .= local_to_owner(mdof_ids)
  for (gid, (lid, owner)) in new_gids
    lid_to_gid[lid] = gid
    lid_to_owner[lid] = owner
  end
  n_global = global_length(mdof_ids)
  return LocalIndices(n_global, rank, lid_to_gid, lid_to_owner)
end

function generate_aggregated_gids(cell_gids, cell_to_fdofs, fdof_is_agg)
  # Create pos/neg local numberings
  fdof_to_posneg, nlpos, nlneg = map(fdof_is_agg) do fdof_is_agg
    npos = 0
    nneg = 0
    fdof_to_posneg = zeros(Int,length(fdof_is_agg))
    for (fdof,is_agg) in enumerate(fdof_is_agg)
      if !is_agg
        npos += 1
        fdof_to_posneg[fdof] = npos
      else
        nneg += 1
        fdof_to_posneg[fdof] = -nneg
      end
    end
    fdof_to_posneg, npos, nneg
  end |> tuple_of_arrays

  # Reindex cell ids
  cell_to_lposneg = map(cell_to_fdofs,fdof_to_posneg) do cell_to_fdofs,fdof_to_posneg
    lazy_map(Broadcasting(Reindex(fdof_to_posneg)),cell_to_fdofs)
  end

  # Generate global master and aggregated dof ids
  mdof_gids, aggdof_gids = GridapDistributed.generate_posneg_gids(
    cell_gids, cell_to_lposneg, nlpos, nlneg
  )

  mdof_to_fdof, aggdof_to_fdof = map(fdof_to_posneg, nlpos, nlneg) do fdof_to_posneg, npos, nneg
    mdof_to_fdof = zeros(Int32,npos)
    aggdof_to_fdof = zeros(Int32,nneg)
    for (fdof,posneg) in enumerate(fdof_to_posneg)
      if posneg > 0
        mdof_to_fdof[posneg] = Int32(fdof)
      elseif posneg < 0
        aggdof_to_fdof[-posneg] = Int32(fdof)
      end
    end
    mdof_to_fdof, aggdof_to_fdof
  end |> tuple_of_arrays

  return mdof_gids, aggdof_gids, mdof_to_fdof, aggdof_to_fdof, fdof_to_posneg
end

function get_aggdof_ptrs(
  aggdof_to_fdof,
  fdof_to_cell,
  cell_to_root,
  cell_to_fdofs,
  aggdof_ids = eachindex(aggdof_to_fdof)
)
  ptrs = zeros(Int32, length(aggdof_to_fdof)+1)
  cache = array_cache(cell_to_fdofs)
  for aggdof in aggdof_ids
    fdof = aggdof_to_fdof[aggdof]
    cell = fdof_to_cell[fdof]
    root = cell_to_root[cell]
    dofs = getindex!(cache,cell_to_fdofs,root)
    ptrs[aggdof+1] = length(dofs)
  end
  return ptrs
end

function aggdof_to_dofs!(
  data, ptrs,
  aggdof_to_fdof,
  fdof_to_cell,
  cell_to_root,
  cell_to_fdofs,
  aggdof_ids = eachindex(aggdof_to_fdof)
)
  cache = array_cache(cell_to_fdofs)
  println(ptrs)
  for aggdof in aggdof_ids
    fdof = aggdof_to_fdof[aggdof]
    cell = fdof_to_cell[fdof]
    root = cell_to_root[cell]
    dofs = getindex!(cache,cell_to_fdofs,root)
    p = ptrs[aggdof]-1
    for (i,dof) in enumerate(dofs)
      data[p+i] = dof
    end
  end
end

function aggdof_to_coeffs!(
  data, ptrs,
  aggdof_to_fdof,
  fdof_to_cell,
  fdof_to_ldof,
  cell_to_coeffs,
  cell_to_proj,
  aggdof_ids = eachindex(aggdof_to_fdof)
)
  z = zero(eltype(eltype(cell_to_coeffs)))
  cache2 = array_cache(cell_to_coeffs)
  cache3 = array_cache(cell_to_proj)

  for aggdof in aggdof_ids
    fdof = aggdof_to_fdof[aggdof]
    cell = fdof_to_cell[fdof]
    coeffs = getindex!(cache2,cell_to_coeffs,cell)
    proj = getindex!(cache3,cell_to_proj,cell)
    ldof = fdof_to_ldof[fdof]
    p = ptrs[aggdof]-1
    for b in axes(proj,2)
      coeff = z
      for c in axes(coeffs,2)
        coeff += coeffs[ldof,c]*proj[c,b]
      end
      data[p+b] = coeff
    end
  end
end
