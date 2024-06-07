
struct DistributedEmbeddedDiscretization{A,B} <: GridapType
  discretizations::A
  model::B
  function DistributedEmbeddedDiscretization(
    d::AbstractArray{<:AbstractEmbeddedDiscretization},
    model::DistributedDiscreteModel)
    A = typeof(d)
    B = typeof(model)
    new{A,B}(d,model)
  end
end

local_views(a::DistributedEmbeddedDiscretization) = a.discretizations

get_background_model(a::DistributedEmbeddedDiscretization) = a.model

function get_geometry(a::DistributedEmbeddedDiscretization)
  loc_geometries = map(get_geometry,local_views(a))
  get_geometry(loc_geometries)
end

function get_geometry(a::AbstractArray{<:CSG.Geometry})
  PartitionedArrays.getany(a)
end

function cut(bgmodel::DistributedDiscreteModel,args...)
  cut(LevelSetCutter(),bgmodel,args...)
end

function cut(cutter::Cutter,bgmodel::DistributedDiscreteModel,args...)
  gids = get_cell_gids(bgmodel)
  cuts = map(local_views(bgmodel),local_views(gids)) do bgmodel,gids
    ownmodel = remove_ghost_cells(bgmodel,gids)
    cutgeo = cut(cutter,ownmodel,args...)
    change_bgmodel(cutgeo,bgmodel,own_to_local(gids))
  end
  consistent_bgcell_to_inoutcut!(cuts,gids)
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

function cut_facets(bgmodel::DistributedDiscreteModel,args...)
  cut_facets(LevelSetCutter(),bgmodel,args...)
end

function cut_facets(cutter::Cutter,bgmodel::DistributedDiscreteModel,args...)
  D = map(num_dims,local_views(bgmodel)) |> PartitionedArrays.getany
  cell_gids = get_cell_gids(bgmodel)
  facet_gids = get_face_gids(bgmodel,D-1)
  cuts = map(
    local_views(bgmodel),
    local_views(cell_gids),
    local_views(facet_gids)) do bgmodel,cell_gids,facet_gids
    ownmodel = remove_ghost_cells(bgmodel,cell_gids)
    facet_to_pfacet = get_face_to_parent_face(ownmodel,D-1)
    cutfacets = cut_facets(cutter,ownmodel,args...)
    cutfacets = change_bgmodel(cutfacets,bgmodel,facet_to_pfacet)
    remove_ghost_subfacets(cutfacets,facet_gids)
  end
  consistent_bgfacet_to_inoutcut!(cuts,facet_gids)
  DistributedEmbeddedDiscretization(cuts,bgmodel)
end

function Triangulation(
  cutgeo::DistributedEmbeddedDiscretization,
  in_or_out::ActiveInOrOut,
  args...)

  distributed_embedded_triangulation(Triangulation,cutgeo,in_or_out,args...)
end

function Triangulation(cutgeo::DistributedEmbeddedDiscretization,args...)
  trian = distributed_embedded_triangulation(Triangulation,cutgeo,args...)
  remove_ghost_cells(trian)
end

function EmbeddedBoundary(cutgeo::DistributedEmbeddedDiscretization,args...)
  trian = distributed_embedded_triangulation(EmbeddedBoundary,cutgeo,args...)
  remove_ghost_cells(trian)
end

function SkeletonTriangulation(cutgeo::DistributedEmbeddedDiscretization,args...)
  trian = distributed_embedded_triangulation(SkeletonTriangulation,cutgeo,args...)
  remove_ghost_cells(trian)
end

function BoundaryTriangulation(cutgeo::DistributedEmbeddedDiscretization,args...)
  trian = distributed_embedded_triangulation(BoundaryTriangulation,cutgeo,args...)
  remove_ghost_cells(trian)
end

function distributed_embedded_triangulation(
  T,
  cutgeo::DistributedEmbeddedDiscretization,
  args...)

  trians = map(local_views(cutgeo)) do lcutgeo
    T(lcutgeo,args...)
  end
  bgmodel = get_background_model(cutgeo)
  DistributedTriangulation(trians,bgmodel)
end

function compute_bgfacet_to_inoutcut(
  bgmodel::DistributedDiscreteModel,
  bgf_to_ioc::AbstractArray{<:AbstractVector})

  D = num_dims(eltype(local_views(bgmodel)))
  gids = get_cell_gids(bgmodel)
  bgf_to_ioc = map(
    local_views(bgmodel),
    local_views(gids),
    bgf_to_ioc) do bgmodel,gids,bgf_to_ioc

    ownmodel = remove_ghost_cells(bgmodel,gids)
    f_to_pf = Gridap.Geometry.get_face_to_parent_face(ownmodel,D-1)
    _bgf_to_ioc = Vector{eltype(bgf_to_ioc)}(undef,num_faces(bgmodel,D-1))
    _bgf_to_ioc[f_to_pf] .= bgf_to_ioc
    _bgf_to_ioc
  end
  facet_gids = get_face_gids(bgmodel,D-1)
  pbgf_to_ioc = PVector(bgf_to_ioc,partition(facet_gids))
  consistent!(pbgf_to_ioc) |> wait
  local_values(pbgf_to_ioc)
end

function compute_bgfacet_to_inoutcut(
  cutter::Cutter,
  bgmodel::DistributedDiscreteModel,
  geo)

  gids = get_cell_gids(bgmodel)
  bgf_to_ioc = map(local_views(bgmodel),local_views(gids)) do model,gids
    ownmodel = remove_ghost_cells(model,gids)
    compute_bgfacet_to_inoutcut(cutter,ownmodel,geo)
  end
  compute_bgfacet_to_inoutcut(bgmodel,bgf_to_ioc)
end

function compute_bgfacet_to_inoutcut(bgmodel::DistributedDiscreteModel,args...)
  cutter = LevelSetCutter()
  compute_bgfacet_to_inoutcut(cutter,bgmodel,args...)
end

function compute_bgcell_to_inoutcut(cutgeo::DistributedEmbeddedDiscretization,args...)
  map(local_views(cutgeo)) do cutgeo
    compute_bgcell_to_inoutcut(cutgeo,args...)
  end
end

function remove_ghost_cells(trian::DistributedTriangulation)
  model = get_background_model(trian)
  gids = get_cell_gids(model)
  trians = map(local_views(trian),local_views(gids)) do trian,gids
    remove_ghost_cells(trian,gids)
  end
  DistributedTriangulation(trians,model)
end

function remove_ghost_cells(trian::AppendedTriangulation,gids)
  a = remove_ghost_cells(trian.a,gids)
  b = remove_ghost_cells(trian.b,gids)
  lazy_append(a,b)
end

function remove_ghost_cells(trian::SubFacetTriangulation,gids)
  model = get_background_model(trian)
  D     = num_cell_dims(model)
  glue  = get_glue(trian,Val{D}())
  remove_ghost_cells(glue,trian,gids)
end

function remove_ghost_cells(model::DiscreteModel,gids::AbstractLocalIndices)
  DiscreteModelPortion(model,own_to_local(gids))
end

function consistent_bgcell_to_inoutcut!(
  cuts::AbstractArray{<:AbstractEmbeddedDiscretization},
  gids::PRange)

  ls_to_bgcell_to_inoutcut = map(get_ls_to_bgcell_to_inoutcut,cuts)
  _consistent!(ls_to_bgcell_to_inoutcut,gids)
end

function get_ls_to_bgcell_to_inoutcut(cut::EmbeddedDiscretization)
  cut.ls_to_bgcell_to_inoutcut
end

function consistent_bgfacet_to_inoutcut!(
  cuts::AbstractArray{<:AbstractEmbeddedDiscretization},
  gids::PRange)

  ls_to_bgfacet_to_inoutcut = map(get_ls_to_bgfacet_to_inoutcut,cuts)
  _consistent!(ls_to_bgfacet_to_inoutcut,gids)
end

function get_ls_to_bgfacet_to_inoutcut(cut::EmbeddedFacetDiscretization)
  cut.ls_to_facet_to_inoutcut
end

function _consistent!(
  p_to_i_to_a::AbstractArray{<:Vector{<:Vector}},
  prange::PRange)

  n = map(length,p_to_i_to_a) |> PartitionedArrays.getany
  for i in 1:n
    p_to_a = map(i_to_a->i_to_a[i],p_to_i_to_a)
    PVector(p_to_a,partition(prange)) |> consistent! |> wait
    map(p_to_a,p_to_i_to_a) do p_to_a,p_to_ia
      copyto!(p_to_ia[i],p_to_a)
    end
  end
end

function change_bgmodel(
  cutgeo::DistributedEmbeddedDiscretization,
  model::DistributedDiscreteModel,
  args...)

  cuts = _change_bgmodels(cutgeo,model,args...)
  gids = get_cell_gids(model)
  ls_to_bgcell_to_inoutcut = map(c->c.ls_to_bgcell_to_inoutcut,cuts)
  _consistent!(ls_to_bgcell_to_inoutcut,gids)
  DistributedEmbeddedDiscretization(cuts,model)
end

function _change_bgmodels(
  cutgeo::DistributedEmbeddedDiscretization,
  model::DistributedDiscreteModel,
  cell_to_newcell)

  map(local_views(cutgeo),local_views(model),cell_to_newcell) do c,m,c_to_nc
    change_bgmodel(c,m,c_to_nc)
  end
end

function _change_bgmodels(
  cutgeo::DistributedEmbeddedDiscretization,
  model::DistributedDiscreteModel)

  map(local_views(cutgeo),local_views(model)) do c,m
    change_bgmodel(c,m)
  end
end

function change_bgmodel(
  cut::EmbeddedDiscretization,
  newmodel::DiscreteModel,
  cell_to_newcell=1:num_cells(get_background_model(cut)))

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
    cut.geo)
end

function change_bgmodel(
  cut::EmbeddedFacetDiscretization,
  newmodel::DiscreteModel,
  facet_to_newfacet=1:num_facets(get_background_model(cut)))

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
    cut.geo)
end

function change_bgmodel(cells::SubCellData,cell_to_newcell)
  cell_to_bgcell = lazy_map(Reindex(cell_to_newcell),cells.cell_to_bgcell)
  SubCellData(
    cells.cell_to_points,
    collect(Int32,cell_to_bgcell),
    cells.point_to_coords,
    cells.point_to_rcoords)
end

function change_bgmodel(facets::SubFacetData,cell_to_newcell)
  facet_to_bgcell = lazy_map(Reindex(cell_to_newcell),facets.facet_to_bgcell)
  SubFacetData(
    facets.facet_to_points,
    facets.facet_to_normal,
    collect(Int32,facet_to_bgcell),
    facets.point_to_coords,
    facets.point_to_rcoords)
end

function remove_ghost_subfacets(cut::EmbeddedFacetDiscretization,facet_gids)
  bgfacet_mask = map(!iszero,local_to_owner(facet_gids))
  subfacet_mask = map(Reindex(bgfacet_mask),cut.subfacets.cell_to_bgcell)
  new_subfacets = findall(subfacet_mask)
  subfacets = SubCellData(cut.subfacets,new_subfacets)
  ls_to_subfacet_to_inout = map(cut.ls_to_subfacet_to_inout) do sf_to_io
    map(Reindex(sf_to_io),new_subfacets)
  end
  EmbeddedFacetDiscretization(
    cut.bgmodel,
    cut.ls_to_facet_to_inoutcut,
    subfacets,
    ls_to_subfacet_to_inout,
    cut.oid_to_ls,
    cut.geo)
end
