
function CellData.Measure(
  t::DistributedTriangulation,
  quad::Tuple{MomentFitted,Vararg};
  kwargs...)
  @notimplemented
  name, _args, _kwargs = quad
  cut,cutfacets,_args... = _args
  t = remove_ghost_cells(t)
  measures = map(
    local_views(t),
    local_views(cut),
    local_views(cutfacets)) do trian,cut,cutfacets
      quad = name, (cut,cutfacets,_args...), _kwargs
      Measure(trian,quad;kwargs...)
  end
  DistributedMeasure(measures)
end
