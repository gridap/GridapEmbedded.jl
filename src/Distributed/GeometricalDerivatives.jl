

function remove_ghost_cells(
  trian::Union{<:CutFaceBoundaryTriangulation,<:CutFaceSkeletonTriangulation},gids
)
  model = get_background_model(trian)
  Dm    = num_cell_dims(model)
  glue  = get_glue(trian,Val(Dm))
  remove_ghost_cells(glue,trian,gids)
end

for func in (:get_subfacet_normal_vector,:get_ghost_normal_vector,:get_conormal_vector)
  @eval begin
    function $func(a::DistributedTriangulation)
      fields = map($func,local_views(a))
      DistributedCellField(fields,a)
    end
  end
end

function LevelSetCutters.DifferentiableTriangulation(trian::DistributedTriangulation,fe_space)
  model = get_background_model(trian)
  trians = map(DifferentiableTriangulation,local_views(trian),local_views(fe_space))
  return DistributedTriangulation(trians,model)
end

function FESpaces._change_argument(
  op,f,
  local_trians::AbstractArray{<:Union{<:DifferentiableTriangulation,<:DifferentiableAppendedTriangulation,<:DifferentiableTriangulationView}},
  uh::GridapDistributed.DistributedADTypes
)
  function dist_cf(uh::DistributedCellField,cfs)
    DistributedCellField(cfs,get_triangulation(uh))
  end
  function dist_cf(uh::DistributedMultiFieldCellField,cfs)
    sf_cfs = map(DistributedCellField,
      [tuple_of_arrays(map(cf -> Tuple(cf.single_fields),cfs))...],
      map(get_triangulation,uh)
    )
    DistributedMultiFieldCellField(sf_cfs,cfs)
  end

  uhs = local_views(uh)
  spaces = map(get_fe_space,uhs)
  function g(cell_u)
    cfs = map(CellField,spaces,cell_u)
    cf = dist_cf(uh,cfs)
    map(update_trian!,local_trians,spaces,local_views(cf))
    cg = f(cf)
    map(local_trians,spaces) do Ω, V
      update_trian!(Ω,V,nothing)
    end
    map(get_contribution,local_views(cg),local_trians)
  end
  g
end
