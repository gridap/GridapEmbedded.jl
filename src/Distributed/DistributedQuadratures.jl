
function CellData.Measure(
  trian::DistributedTriangulation,
  quad::Tuple{MomentFitted,Vararg};
  kwargs...)
  @notimplemented
  name, _args, _kwargs = quad
  cut,cutfacets,_args... = _args

  model = get_background_model(trian)
  gids  = get_cell_gids(model)
  trian = remove_ghost_cells(trian,gids)
  measures = map(local_views(trian),local_views(cut),local_views(cutfacets)) do trian,cut,cutfacets
    quad = name, (cut,cutfacets,_args...), _kwargs
    Measure(trian,quad;kwargs...)
  end
  DistributedMeasure(measures)
end
