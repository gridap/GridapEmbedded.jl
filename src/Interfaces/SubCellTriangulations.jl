
struct SubCellData{Dr,Dp,Tp,Tr} <: GridapType
  cell_to_points::Table{Int32,Vector{Int32},Vector{Int32}}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,Tp}}
  point_to_rcoords::Vector{Point{Dr,Tr}}
end

function SubCellData(st::SubCellData,newcells::AbstractVector{<:Integer})
  cell_to_points = Table(lazy_map(Reindex(st.cell_to_points),newcells))
  cell_to_bgcell = st.cell_to_bgcell[newcells]
  SubCellData(
    cell_to_points,
    cell_to_bgcell,
    st.point_to_coords,
    st.point_to_rcoords)
end

# Implementation of Triangulation interface

struct SubCellTriangulation{Dc,Dp,T,A} <: Triangulation{Dc,Dp}
  subcells::SubCellData{Dc,Dp,T}
  bgmodel::A
  subgrid::UnstructuredGrid{Dc,Dp,T,NonOriented,Nothing}
  function SubCellTriangulation(
    subcells::SubCellData{Dc,Dp,T},bgmodel::DiscreteModel) where {Dc,Dp,T}
    subgrid = UnstructuredGrid(subcells)
    A = typeof(bgmodel)
    new{Dc,Dp,T,A}(subcells,bgmodel,subgrid)
  end
end

function SubCellTriangulation(subcells::SubCellData,bgmodel::AdaptedDiscreteModel) 
  SubCellTriangulation(subcells,bgmodel.model)
end

function get_background_model(a::SubCellTriangulation)
  a.bgmodel
end

function get_active_model(a::SubCellTriangulation)
  msg = """
  This is not implemented, but also not needed in practice.
  Embedded Grids implemented for integration, not interpolation.
  """
  @notimplemented  msg
end

function get_grid(a::SubCellTriangulation)
  a.subgrid
end

function get_glue(a::SubCellTriangulation{Dc},::Val{D}) where {Dc,D}
  if D != Dc
    return nothing
  end
  tface_to_mface = a.subcells.cell_to_bgcell
  tface_to_mface_map = _setup_cell_ref_map(a.subcells,a.subgrid)
  FaceToFaceGlue(tface_to_mface,tface_to_mface_map,nothing)
end

function move_contributions(scell_to_val::AbstractArray,strian::SubCellTriangulation)
  model = get_background_model(strian)
  ncells = num_cells(model)
  cell_to_touched = fill(false,ncells)
  scell_to_cell = strian.subcells.cell_to_bgcell
  cell_to_touched[scell_to_cell] .= true
  Ωa = Triangulation(model,cell_to_touched)
  acell_to_val = move_contributions(scell_to_val,strian,Ωa)
  acell_to_val, Ωa 
end

function _setup_cell_ref_map(st,grid)
  cell_to_points = st.cell_to_points
  point_to_rcoords = st.point_to_rcoords
  cell_to_rcoords = lazy_map(Broadcasting(Reindex(point_to_rcoords)),cell_to_points)
  ctype_to_reffe = get_reffes(grid)
  cell_to_ctype = get_cell_type(grid)
  @notimplementedif length(ctype_to_reffe) != 1
  reffe = first(ctype_to_reffe)
  cell_to_shapefuns = expand_cell_data([get_shapefuns(reffe)],cell_to_ctype)
  cell_to_ref_map = lazy_map(linear_combination,cell_to_rcoords,cell_to_shapefuns)
  cell_to_ref_map
end

# API

function UnstructuredGrid(st::SubCellData{D}) where D
  reffe = LagrangianRefFE(Float64,Simplex(Val{D}()),1)
  cell_types = fill(Int8(1),length(st.cell_to_points))
  UnstructuredGrid(
    st.point_to_coords,
    st.cell_to_points,
    [reffe,],
    cell_types)
end

function Visualization.visualization_data(st::SubCellData,filename::String;celldata=Dict())
  ug = UnstructuredGrid(st)
  degree = 0
  trian = GenericTriangulation(ug)
  quad = CellQuadrature(trian,degree)
  dV = integrate(1,quad)
  newcelldata = ["bgcell"=>st.cell_to_bgcell, "dV"=>dV]
  _celldata = Dict()
  for (k,v) in celldata
    _celldata[k] = v
  end
  for (k,v) in newcelldata
    _celldata[k] = v
  end
  (Visualization.VisualizationData(ug,filename,celldata=_celldata),)
end
