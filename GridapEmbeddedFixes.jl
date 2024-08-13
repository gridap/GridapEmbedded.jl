function Gridap.Geometry.get_background_model(
    a::GridapEmbedded.Interfaces.SubCellTriangulation{Dc,
                                                      Dp,
                                                      T,
                                                      <:Gridap.Adaptivity.AdaptedDiscreteModel}) where {Dc,Dp,T}
 Gridap.Adaptivity.get_model(a.bgmodel)
end

const SFC_ADM{Dc,Dp,T} =
    GridapEmbedded.Interfaces.SubFacetTriangulation{Dc,Dp,T,<:Gridap.Adaptivity.AdaptedDiscreteModel}

function Gridap.Geometry.get_background_model(
 a::SFC_ADM{Dc,Dp,T,<:Gridap.Adaptivity.AdaptedDiscreteModel}) where {Dc,Dp,T}
 Gridap.Adaptivity.get_model(a.bgmodel)
end

const APP_ST_AT{Dc,Dp} =  Gridap.Geometry.AppendedTriangulation{Dc,Dp,
                                                               <:GridapEmbedded.Interfaces.SubCellTriangulation,
                                                               <:Gridap.Adaptivity.AdaptedTriangulation}

const APP_AT_AT{Dc,Dp} =  Gridap.Geometry.AppendedTriangulation{Dc,Dp,
                                                               <:Gridap.Adaptivity.AdaptedTriangulation,
                                                               <:Gridap.Adaptivity.AdaptedTriangulation}

function Gridap.FESpaces.get_cell_fe_data(fun,f,ttrian::Union{APP_ST_AT{Dc,Dp},APP_AT_AT{Dc,Dp}}) where {Dc,Dp}
  sface_to_data = fun(f)
  strian = get_triangulation(f)
  Gridap._get_cell_fe_data(fun,sface_to_data, strian, ttrian)
end

function Gridap.FESpaces.get_cell_fe_data(fun,f,ttrian::SFC_ADM{Dc,Dp,T}) where {Dc,Dp,T}
  sface_to_data = fun(f)
  strian = get_triangulation(f)
  Gridap._get_cell_fe_data(fun,sface_to_data, strian, ttrian)
end

function  Gridap._get_cell_fe_data(fun,
  sface_to_data,
  strian::Gridap.Geometry.BodyFittedTriangulation{Dc,Dp},
  ttrian::Union{APP_ST_AT{Dc,Dp},APP_AT_AT{Dc,Dp}}) where {Dc,Dp}

  bt_adapted_model = Gridap.Adaptivity.get_adapted_model(ttrian.b)
  if (get_background_model(strian) === Gridap.Adaptivity.get_parent(bt_adapted_model))
    tglue = get_glue(ttrian,Val(Dc))
    ttrian_tface_to_mface=tglue.tface_to_mface
    strian_ttrian_tface_to_mface=bt_adapted_model.glue.n2o_faces_map[end][ttrian_tface_to_mface]
    sglue = get_glue(strian,Val(Dc))
    strian_mface_to_tface=sglue.mface_to_tface
    lazy_map(Reindex(sface_to_data),strian_mface_to_tface[strian_ttrian_tface_to_mface])
  elseif (get_background_model(strian) === get_background_model(bt))
    Gridap.Helpers.@notimplemented
  else
    Gridap.Helpers.@notimplemented
  end
end

function Gridap._get_cell_fe_data(fun,
  sface_to_data,
  strian::Gridap.Geometry.BodyFittedTriangulation{Dc,Dp},
  ttrian::SFC_ADM{Df,Dp}) where {Dc,Df,Dp}

  if (get_background_model(strian) === Gridap.Adaptivity.get_parent(ttrian.bgmodel))
    tcell_to_mcell=ttrian.subfacets.facet_to_bgcell
    adapt_glue=ttrian.bgmodel.glue
    s_coarse_mfaces=adapt_glue.n2o_faces_map[end][tcell_to_mcell]
    sglue=get_glue(strian,Val(Dc))
    s_mface_to_tface=sglue.mface_to_tface
    sface_to_data_reindexed=lazy_map(Reindex(sface_to_data),s_mface_to_tface[s_coarse_mfaces])
  else
    Gridap.Helpers.@notimplemented
  end
end

function Gridap.Geometry.move_contributions(scell_to_val::AbstractArray,
                                            strian::GridapEmbedded.Interfaces.SubFacetTriangulation)
  model = strian.bgmodel # We cannot call get_background_model(strian) because it does not return the adapted model!!!
  ncells = num_cells(model)
  cell_to_touched = fill(false,ncells)
  scell_to_cell = strian.subfacets.facet_to_bgcell
  cell_to_touched[scell_to_cell] .= true
  Ωa = Triangulation(model,cell_to_touched)
  acell_to_val = move_contributions(scell_to_val,strian,Ωa)
  acell_to_val, Ωa
end

function Gridap.Geometry.move_contributions(scell_to_val::AbstractArray,
                                            strian::GridapEmbedded.Interfaces.SubCellTriangulation)
  model = strian.bgmodel # We cannot call get_background_model(strian) because it does not return the adapted model!!!
  ncells = num_cells(model)
  cell_to_touched = fill(false,ncells)
  scell_to_cell = strian.subcells.cell_to_bgcell
  cell_to_touched[scell_to_cell] .= true
  Ωa = Triangulation(model,cell_to_touched)
  acell_to_val = move_contributions(scell_to_val,strian,Ωa)
  acell_to_val, Ωa
end

"""
- strian: denoted by ΩHact, the coarse active triangulation for the FE space
- ttrian: denoted by Ωh, the refined physical triangulation for integration

ΩHact is created by deactivating some cells on the coarse mesh mH.
Ωh consists of cut cells (Ωh.a) and interior cells (Ωh.b), both are parts of the refined mesh mh.
The cell ids mapping required is from Ωh to ΩHact.
We can achieve this by the path: Ωh -> mh -> mH -> ΩHact
"""
function Gridap.FESpaces._compute_cell_ids(uh, strian::Triangulation, ttrian::APP_ST_AT)
  bbgmodel = Gridap.Adaptivity.get_adapted_model(ttrian.b)
  Gridap.Helpers.@check get_background_model(strian) === Gridap.Adaptivity.get_parent(bbgmodel)
  D = num_cell_dims(strian)
  sglue, tglue = get_glue(strian, Val(D)), get_glue(ttrian, Val(D))
  cell_id = Vector{Int32}(undef, num_cells(ttrian))
  ωh_to_mh = tglue.tface_to_mface              # Ωh -> mh
  mh_to_mH = bbgmodel.glue.n2o_faces_map[end]  # mh -> mH
  ωH_to_mH = sglue.mface_to_tface              # mH -> ΩHact
  for (ia, ic) in enumerate(ωh_to_mh)
    cell_id[ia] = ωH_to_mH[mh_to_mH[ic]]
  end
  cell_id
end
