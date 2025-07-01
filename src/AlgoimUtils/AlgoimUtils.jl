module AlgoimUtils

using CxxWrap
using Algoim
using Gridap

import Base.prod
@inline prod(x::Point{N,T}) where {N,T} = prod(x.data)

import Algoim: to_array
@inline to_array(x::Point{N,T}) where {N,T} = collect(x.data)

using Gridap.Helpers
using Gridap.ReferenceFEs
using Gridap.Arrays: collect1d, CompressedArray, Table
using Gridap.Adaptivity
using Gridap.Fields: testitem
using Gridap.CellData
using Gridap.CellData: GenericField, GenericCellField, get_data
using Gridap.CellData: _point_to_cell_cache, _point_to_cell!
using Gridap.Geometry
using Gridap.FESpaces

using PartitionedArrays
using GridapDistributed
using GridapDistributed: DistributedDiscreteModel
using GridapDistributed: DistributedCartesianDiscreteModel
using GridapDistributed: DistributedAdaptedDiscreteModel
using GridapDistributed: DistributedTriangulation
using GridapDistributed: DistributedMeasure
using GridapDistributed: DistributedCellField
using GridapDistributed: DistributedFESpace
using GridapDistributed: DistributedSingleFieldFEFunction
using GridapDistributed: DistributedVisualizationData
using GridapDistributed: add_ghost_cells
import GridapDistributed: local_views

using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: Simplex

using MiniQhull
using FillArrays

import Algoim: CachedLevelSetValue
import Algoim: CachedLevelSetGradient
import Algoim: AlgoimCallLevelSetFunction
import Algoim: normal
import Gridap.ReferenceFEs: Quadrature

export CellQuadratureAndActiveMask
export TriangulationAndMeasure
export algoim
export Quadrature
export is_cell_active
export restrict_measure
export fill_cpp_data
export compute_closest_point_projections
export compute_normal_displacement
export compute_normal_displacement!
export compute_distance_fe_function
export active_triangulation
export narrow_band_triangulation
export delaunaytrian
export convexhull

struct Algoim <: QuadratureName end
const algoim = Algoim()

(φ::CachedLevelSetValue{<:CellField})(p,i::Float32) = begin
  _p = Point(p)
  (carr,cgetindex,ceval) = φ.c
  evaluate!(ceval,getindex!(cgetindex,carr,Int(i)),_p)
end

(φ::CachedLevelSetGradient{<:CellField})(p,i::Float32) = begin
  _p = Point(p)
  (carr,cgetindex,ceval) = φ.c
  _val = evaluate!(ceval,getindex!(cgetindex,carr,Int(i)),_p)
  ConstCxxRef(to_uvector(to_array(_val)))
end

function AlgoimCallLevelSetFunction(φ::CellField,∇φ::CellField)
  φtrian = get_triangulation(φ)
  ∇φtrian = get_triangulation(∇φ)
  xφ = testitem(testitem(get_cell_coordinates(φtrian)))
  x∇φ = testitem(testitem(get_cell_coordinates(∇φtrian)))
  cellφ = get_array(φ)
  cell∇φ = get_array(∇φ)
  cache_φ = (cellφ,array_cache(cellφ),return_cache(testitem(cellφ),xφ))
  cache_∇φ = (cell∇φ,array_cache(cell∇φ),return_cache(testitem(cell∇φ),x∇φ))
  AlgoimCallLevelSetFunction{typeof(φ),typeof(∇φ),typeof(cache_φ),typeof(cache_∇φ)}(φ,∇φ,cache_φ,cache_∇φ)
end

struct DistributedAlgoimCallLevelSetFunction{A<:AbstractArray{<:AlgoimCallLevelSetFunction}} <: GridapType
  values   ::DistributedCellField
  gradients::DistributedCellField
  levelsets::A
end

local_views(a::DistributedAlgoimCallLevelSetFunction) = a.levelsets

function AlgoimCallLevelSetFunction(φ::DistributedCellField,∇φ::DistributedCellField)
  levelsets = map((v,g)->AlgoimCallLevelSetFunction(v,g),local_views(φ),local_views(∇φ))
  DistributedAlgoimCallLevelSetFunction(φ,∇φ,levelsets)
end

function normal(phi::AlgoimCallLevelSetFunction,x::AbstractVector{<:Point},cell_id::Int=1)
  map(xi->normal(phi,xi,cell_id),x)
end

function normal(phi::AlgoimCallLevelSetFunction,trian::Triangulation)
  f = [ x -> normal(phi,x,i) for i in 1:num_cells(trian) ]
  f = map(GenericField,f)
  GenericCellField(f,trian,PhysicalDomain())
end

function normal(phi::AlgoimCallLevelSetFunction,trian::DistributedTriangulation)
  normals = map(t->normal(phi,t),local_views(trian))
  DistributedCellField(normals,trian)
end

function normal(phi::DistributedAlgoimCallLevelSetFunction,_trian::DistributedTriangulation)
  trian = add_ghost_cells(_trian)
  normals = map((φ,t)->normal(φ,t),local_views(phi),local_views(trian))
  DistributedCellField(normals,trian)
end

function normal(ls::AlgoimCallLevelSetFunction{<:CellField,<:CellField},x::Point,cell_id::Int=1)
  (carr,cgetindex,ceval) = ls.cache_∇φ
  gx = evaluate!(ceval,getindex!(cgetindex,carr,cell_id),x)
  gx/norm(gx)
end

function Quadrature(trian::Grid,::Algoim,phi::LevelSetFunction,degree::Int;kwargs...)
  ctype_polytope = map(get_polytope,get_reffes(trian))
  @notimplementedif !all(map(is_n_cube,ctype_polytope))
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_bboxes = collect1d(lazy_map(a->(a[1],a[end]),cell_to_coords))
  jls = JuliaFunctionLevelSet(phi,Val{num_dims(trian)}())
  cell_to_quad = map(enumerate(cell_to_bboxes)) do (cell_id,bbox)
    bbmin, bbmax = bbox
    Quadrature(cell_id,bbmin,bbmax,jls,phi,degree;kwargs...)
  end
  CompressedArray(cell_to_quad,1:length(cell_to_quad))
end

function Quadrature(trian::Grid,::Algoim,phi::LevelSetFunction,
                    own_to_local::AbstractVector,degree::Int;kwargs...)
  ctype_polytope = map(get_polytope,get_reffes(trian))
  @notimplementedif !all(map(is_n_cube,ctype_polytope))
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_bboxes = collect1d(lazy_map(a->(a[1],a[end]),cell_to_coords))
  jls = JuliaFunctionLevelSet(phi,Val{num_dims(trian)}())
  cell_to_quad = map(enumerate(cell_to_bboxes)) do (own_cell_id,bbox)
    bbmin, bbmax = bbox
    cell_id = own_to_local[own_cell_id]
    Quadrature(cell_id,bbmin,bbmax,jls,phi,degree;kwargs...)
  end
  CompressedArray(cell_to_quad,1:length(cell_to_quad))
end

function Quadrature(trian::Grid,::Algoim,
                    phi1::LevelSetFunction,phi2::LevelSetFunction,
                    degree::Int;kwargs...)
  ctype_polytope = map(get_polytope,get_reffes(trian))
  @notimplementedif !all(map(is_n_cube,ctype_polytope))
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_bboxes = collect1d(lazy_map(a->(a[1],a[end]),cell_to_coords))
  jls1 = JuliaFunctionLevelSet(phi1,Val{num_dims(trian)}())
  jls2 = JuliaFunctionLevelSet(phi2,Val{num_dims(trian)}())
  cell_to_quad = map(enumerate(cell_to_bboxes)) do (cell_id,bbox)
    bbmin, bbmax = bbox
    Quadrature(cell_id,bbmin,bbmax,jls1,jls2,phi1,phi2,degree;kwargs...)
  end
  CompressedArray(cell_to_quad,1:length(cell_to_quad))
end

function Quadrature(trian::Grid,::Algoim,
                    phi1::LevelSetFunction,phi2::LevelSetFunction,
                    own_to_local::AbstractVector,degree::Int;kwargs...)
  ctype_polytope = map(get_polytope,get_reffes(trian))
  @notimplementedif !all(map(is_n_cube,ctype_polytope))
  cell_to_coords = get_cell_coordinates(trian)
  cell_to_bboxes = collect1d(lazy_map(a->(a[1],a[end]),cell_to_coords))
  jls1 = JuliaFunctionLevelSet(phi1,Val{num_dims(trian)}())
  jls2 = JuliaFunctionLevelSet(phi2,Val{num_dims(trian)}())
  cell_to_quad = map(enumerate(cell_to_bboxes)) do (own_cell_id,bbox)
    bbmin, bbmax = bbox
    cell_id = own_to_local[own_cell_id]
    Quadrature(cell_id,bbmin,bbmax,jls1,jls2,phi1,phi2,degree;kwargs...)
  end
  CompressedArray(cell_to_quad,1:length(cell_to_quad))
end

function Quadrature(cell_id::Int,
                    xmin::Point{N,T},
                    xmax::Point{N,T},
                    jls::LevelSetFunction,
                    phi::LevelSetFunction,
                    degree::Int;
                    phase::Int=CUT) where {N,T}
  coords, weights = fill_quad_data(jls,phi,xmin,xmax,phase,degree,cell_id)
  GenericQuadrature(coords,weights,"Algoim quadrature of degree $degree")
end

function Quadrature(cell_id::Int,
                    xmin::Point{N,T},
                    xmax::Point{N,T},
                    jls1::LevelSetFunction,
                    jls2::LevelSetFunction,
                    phi1::LevelSetFunction,
                    phi2::LevelSetFunction,
                    degree::Int;
                    phase1::Int=IN,
                    phase2::Int=IN) where {N,T}
  coords, weights = fill_quad_data(jls1,jls2,phi1,phi2,xmin,xmax,
                                   phase1,phase2,degree,cell_id)
  GenericQuadrature(coords,weights,"Algoim quadrature of degree $degree")
end

function is_cell_active(cell_quad::AbstractArray{<:Quadrature})
  has_non_empty_quad(x) = num_points(x) > 0
  lazy_map(has_non_empty_quad,cell_quad)
end

function is_cell_active(meas::Measure)
  is_cell_active(get_data(meas.quad))
end

function is_cell_active(meas::DistributedMeasure)
  map(local_views(meas)) do meas
    is_cell_active(get_data(meas.quad))
  end
end

function _cell_quadrature_and_active_mask(trian::Grid,
                                          ::Algoim,phi,args;kwargs)
  cell_quad = Quadrature(trian,algoim,phi,args...;kwargs...)
  cell_to_is_active = is_cell_active(cell_quad)
  cell_quad, cell_to_is_active
end

function _cell_quadrature_and_active_mask(trian::Grid,
                                          ::Algoim,phi1,phi2,args;kwargs)
  cell_quad = Quadrature(trian,algoim,phi1,phi2,args...;kwargs...)
  cell_to_is_active = is_cell_active(cell_quad)
  cell_quad, cell_to_is_active
end

function _cell_quadrature_and_active_mask(trian::DistributedTriangulation,
    ::Algoim,
    phi::AlgoimCallLevelSetFunction,
    args;kwargs)
  ltrians = local_views(trian)
  cell_quad = map(t->Quadrature(t,algoim,phi,args...;kwargs...),ltrians)
  cell_to_is_active = map(cq->is_cell_active(cq),cell_quad)
  cell_quad, cell_to_is_active
end

function _cell_quadrature_and_active_mask(trian::DistributedTriangulation,
    ::Algoim,
    phi1::AlgoimCallLevelSetFunction,
    phi2::AlgoimCallLevelSetFunction,
    args;kwargs)
  ltrians = local_views(trian)
  cell_quad = map(t->Quadrature(t,algoim,phi1,phi2,args...;kwargs...),ltrians)
  cell_to_is_active = map(cq->is_cell_active(cq),cell_quad)
  cell_quad, cell_to_is_active
end

function _cell_quadrature_and_active_mask(trian::DistributedTriangulation,
    ::Algoim,
    phi::DistributedAlgoimCallLevelSetFunction,
    args;kwargs)
  ltrians = local_views(trian); lphis = local_views(phi)
  phitrian = get_triangulation(phi.values)
  gids = get_cell_gids(get_background_model(phitrian))
  own_to_local = map(local_views(phitrian),local_views(gids)) do t,g
    findall(!iszero,local_to_own(g)[t.tface_to_mface])
  end
  cell_quad = map(
    (t,p,otl)->Quadrature(t,algoim,p,otl,args...;kwargs...),ltrians,lphis,own_to_local)
  cell_to_is_active = map(cq->is_cell_active(cq),cell_quad)
  cell_quad, cell_to_is_active
end

function _cell_quadrature_and_active_mask(trian::DistributedTriangulation,
    ::Algoim,
    phi1::DistributedAlgoimCallLevelSetFunction,
    phi2::AlgoimCallLevelSetFunction,
    args;kwargs)
  ltrians = local_views(trian); lphis = local_views(phi1)
  phitrian = get_triangulation(phi1.values)
  gids = get_cell_gids(get_background_model(phitrian))
  own_to_local = map(local_views(phitrian),local_views(gids)) do t,g
    findall(!iszero,local_to_own(g)[t.tface_to_mface])
  end
  cell_quad = map(
    (t,p,otl)->Quadrature(t,algoim,p,phi2,otl,args...;kwargs...),ltrians,lphis,own_to_local)
  cell_to_is_active = map(cq->is_cell_active(cq),cell_quad)
  cell_quad, cell_to_is_active
end

function _cell_quadrature_and_active_mask(trian::DistributedTriangulation,
    ::Algoim,
    phi1::DistributedAlgoimCallLevelSetFunction,
    phi2::DistributedAlgoimCallLevelSetFunction,
    args;kwargs)
  ltrians = local_views(trian); lphis1 = local_views(phi1); lphis2 = local_views(phi2)
  phitrian = get_triangulation(phi1.values)
  gids = get_cell_gids(get_background_model(phitrian))
  own_to_local = map(local_views(phitrian),local_views(gids)) do t,g
    findall(!iszero,local_to_own(g)[t.tface_to_mface])
  end
  cell_quad = map(
    (t,p1,p2,otl)->Quadrature(t,algoim,p1,p2,otl,args...;kwargs...),ltrians,lphis1,lphis2,own_to_local)
  cell_to_is_active = map(cq->is_cell_active(cq),cell_quad)
  cell_quad, cell_to_is_active
end

function _cell_quadrature_and_active_mask(trian,
    name,_args::Tuple{<:AlgoimCallLevelSetFunction,Any};kwargs)
  phi, args = _args
  _cell_quadrature_and_active_mask(trian,name,phi,args;kwargs)
end

function _cell_quadrature_and_active_mask(trian,
    name,_args::Tuple{<:DistributedAlgoimCallLevelSetFunction,Any};kwargs)
  phi, args = _args
  _cell_quadrature_and_active_mask(trian,name,phi,args;kwargs)
end

function _cell_quadrature_and_active_mask(trian,
    name,_args::Tuple{<:DistributedAlgoimCallLevelSetFunction,
                      <:AlgoimCallLevelSetFunction,Any};kwargs)
  phi1, phi2, args = _args
  _cell_quadrature_and_active_mask(trian,name,phi1,phi2,args;kwargs)
end

function _cell_quadrature_and_active_mask(trian,
    name,_args::Tuple{<:DistributedAlgoimCallLevelSetFunction,
                      <:DistributedAlgoimCallLevelSetFunction,Any};kwargs)
  phi1, phi2, args = _args
  _cell_quadrature_and_active_mask(trian,name,phi1,phi2,args;kwargs)
end

function CellQuadratureAndActiveMask(trian,
    quad::Tuple{<:QuadratureName,Any,Any})
  name, args, kwargs = quad
  _cell_quadrature_and_active_mask(trian,name,args;kwargs)
end

function CellQuadratureAndActiveMask(
    model::Union{DiscreteModel,DistributedDiscreteModel},
    quad::Tuple{<:QuadratureName,Any,Any})
  CellQuadratureAndActiveMask(Triangulation(model),quad)
end

function restrict_measure(meas::Measure)
  is_r = is_cell_active(meas)
  rcell_to_cell = findall(is_r)
  oquad = meas.quad
  cell_quad   = lazy_map(Reindex(oquad.cell_quad),rcell_to_cell)
  cell_point  = lazy_map(Reindex(oquad.cell_point),rcell_to_cell)
  cell_weight = lazy_map(Reindex(oquad.cell_weight),rcell_to_cell)
  dds = oquad.data_domain_style
  ids = oquad.integration_domain_style
  trian = Triangulation(oquad.trian,is_r)
  Measure(CellQuadrature(cell_quad,cell_point,cell_weight,trian,dds,ids))
end

function restrict_measure(meas::Measure,trian::Triangulation)
  oquad = meas.quad
  cell_quad   = lazy_map(Reindex(oquad.cell_quad),trian.tface_to_mface)
  cell_point  = lazy_map(Reindex(oquad.cell_point),trian.tface_to_mface)
  cell_weight = lazy_map(Reindex(oquad.cell_weight),trian.tface_to_mface)
  dds = oquad.data_domain_style
  ids = oquad.integration_domain_style
  Measure(CellQuadrature(cell_quad,cell_point,cell_weight,trian,dds,ids))
end

function restrict_measure(meas::Measure,trian::Triangulation,gids::AbstractLocalIndices)
  oquad = meas.quad
  tface_to_own_mface = local_to_own(gids)[trian.tface_to_mface]
  cell_quad   = lazy_map(Reindex(oquad.cell_quad),tface_to_own_mface)
  cell_point  = lazy_map(Reindex(oquad.cell_point),tface_to_own_mface)
  cell_weight = lazy_map(Reindex(oquad.cell_weight),tface_to_own_mface)
  dds = oquad.data_domain_style
  ids = oquad.integration_domain_style
  Measure(CellQuadrature(cell_quad,cell_point,cell_weight,trian,dds,ids))
end

function _triangulation_and_measure(Ωbg::Triangulation,quad::Tuple)
  dΩbg = Measure(Ωbg,quad,data_domain_style=PhysicalDomain())
  # RMK: This is a hack, but algoim interface does not let you 
  # know if a (given) cell intersects the interior of the level 
  # set. The hack consists in inferring this from the size of 
  # each quadrature. I do not expect this hack implies a lot of 
  # extra operations with regards to the proper way to do it.
  cell_to_is_active = is_cell_active(dΩbg)
  Ωᵃ = Triangulation(Ωbg,cell_to_is_active)
  dΩᵃ = restrict_measure(dΩbg,Ωᵃ)
  Ωᵃ,dΩᵃ,cell_to_is_active
end

@inline _measure(Ωbg,quad,phi::AlgoimCallLevelSetFunction) =
  map(lt->Measure(lt,quad,data_domain_style=PhysicalDomain()),local_views(Ωbg))

@inline function _measure(Ωbg,quad,phi::DistributedAlgoimCallLevelSetFunction)
  name, _args, kwargs = quad
  _, args = _args
  ltrians = local_views(Ωbg); lphis = local_views(phi)
  map((lt,lp)->Measure(lt,(name,(lp,args),kwargs),data_domain_style=PhysicalDomain()),ltrians,lphis)
end

function _triangulation_and_measure(Ωbg::Triangulation,
                cell_quad::AbstractArray{<:Quadrature},
                cell_to_is_active::AbstractArray{Bool},
                cell_to_has_non_empty_quad::AbstractArray{Bool})
  dds = PhysicalDomain(); ids = PhysicalDomain()
  dΩbg = Measure(Ωbg,cell_quad,dds,ids)
  Ωᵃ = Triangulation(Ωbg,cell_to_is_active)
  dΩᵃ = restrict_measure(dΩbg,Triangulation(Ωbg,cell_to_has_non_empty_quad))
  Ωᵃ,dΩᵃ,cell_to_is_active
end

function _triangulation_and_measure(Ωbg::DistributedTriangulation,
                           cell_quad::AbstractArray,
                           cell_to_is_active::AbstractArray,
                           cell_to_has_non_empty_quad::AbstractArray)
  ltrians = local_views(Ωbg)
  dds = PhysicalDomain(); ids = PhysicalDomain()
  dΩbg = map((lt,cq)->Measure(lt,cq,dds,ids),ltrians,cell_quad)
  gids = local_views(get_cell_gids(get_background_model(Ωbg)))
  Ωᵃ = map((lt,ca)->Triangulation(lt,ca),ltrians,cell_to_is_active)
  Ωᵖ = map((lt,cp)->Triangulation(lt,cp),ltrians,cell_to_has_non_empty_quad)
  dΩᵃ = map((db,at,gs)->restrict_measure(db,at,gs),dΩbg,Ωᵖ,gids)
  Mbg = get_background_model(Ωbg)
  DΩᵃ = DistributedTriangulation(Ωᵃ,Mbg)
  DΩᵃ,DistributedMeasure(dΩᵃ,DΩᵃ),cell_to_is_active
end

function TriangulationAndMeasure(Ωbg,quad::Tuple)
  msg = "TriangulationAndMeasure can only receive the background triangulation"
  @notimplementedif num_cells(get_background_model(Ωbg)) != num_cells(Ωbg) msg
  _triangulation_and_measure(Ωbg,quad)
end

function TriangulationAndMeasure(Ωbg,cell_quad,cell_to_is_active,
      cell_to_has_non_empty_quad=cell_to_is_active)
  msg = "TriangulationAndMeasure can only receive the background triangulation"
  @notimplementedif num_cells(get_background_model(Ωbg)) != num_cells(Ωbg) msg
  _triangulation_and_measure(Ωbg,
                             cell_quad,
                             cell_to_is_active,
                             cell_to_has_non_empty_quad)
end

using Gridap.Geometry: get_cell_to_parent_cell
using Gridap.CellData: get_cell_quadrature

function compute_closest_point_projections(Ω::Triangulation,φ;
    cppdegree::Int=2,trim::Bool=false,limitstol::Float64=1.0e-8)
  compute_closest_point_projections(get_background_model(Ω),φ,
    cppdegree=cppdegree,trim=trim,limitstol=limitstol)
end

function compute_closest_point_projections(model::CartesianDiscreteModel,
                                           φ::AlgoimCallLevelSetFunction;
                                           cppdegree::Int=2,
                                           trim::Bool=false,
                                           limitstol::Float64=1.0e-8)
  cdesc = get_cartesian_descriptor(model)
  partition = Int32[cdesc.partition...]
  xmin = cdesc.origin
  xmax = xmin + Point(cdesc.sizes .* partition)
  fill_cpp_data(φ,partition,xmin,xmax,cppdegree,trim,limitstol)
end

function compute_closest_point_projections(model::DistributedCartesianDiscreteModel,
                                           φ::AlgoimCallLevelSetFunction;
                                           cppdegree::Int=2,
                                           trim::Bool=false,
                                           limitstol::Float64=1.0e-8)
  gdesc = model.metadata.descriptor
  xmin = gdesc.origin
  gpartition = Int32[gdesc.partition...]
  xmax = xmin + Point(gdesc.sizes .* gpartition)
  cpps = map(local_views(model)) do m
    cdesc = get_cartesian_descriptor(m)
    cmin = ( cdesc.origin - gdesc.origin )
    cmin = round.(cmin.data ./ gdesc.sizes)
    cmin = Int32[cmin...]
    cpartition = Int32[cdesc.partition...]
    cmax = cmin + cpartition
    fill_cpp_data(φ,gpartition,xmin,xmax,
                  cppdegree,trim,limitstol,
                  rmin=cmin,rmax=cmax)
  end
  node_gids = get_face_gids(model,0)
  PVector(cpps,partition(node_gids)) |> consistent! |> wait
  cpps
end

function compute_closest_point_projections(model::DistributedCartesianDiscreteModel,
                                           φ::DistributedAlgoimCallLevelSetFunction;
                                           cppdegree::Int=2,
                                           trim::Bool=false,
                                           limitstol::Float64=1.0e-8)
  gdesc = model.metadata.descriptor
  xmin = gdesc.origin
  gpartition = Int32[gdesc.partition...]
  xmax = xmin + Point(gdesc.sizes .* gpartition)
  cpps = map(local_views(model),local_views(φ)) do m,f
    cdesc = get_cartesian_descriptor(m)
    cmin = ( cdesc.origin - gdesc.origin )
    cmin = round.(cmin.data ./ gdesc.sizes)
    cmin = Int32[cmin...]
    cpartition = Int32[cdesc.partition...]
    cmax = cmin + cpartition
    fill_cpp_data(f,gpartition,xmin,xmax,
                  cppdegree,trim,limitstol,
                  rmin=cmin,rmax=cmax)
  end
  node_gids = get_face_gids(model,0)
  PVector(cpps,partition(node_gids)) |> consistent! |> wait
  cpps
end

function compute_closest_point_projections(fespace::FESpace,
    φ::Union{AlgoimCallLevelSetFunction,DistributedAlgoimCallLevelSetFunction},
    order::Int;cppdegree::Int=2,trim::Bool=false,limitstol::Float64=1.0e-8)
  trian = get_triangulation(fespace)
  model = get_background_model(trian)
  cps = compute_closest_point_projections(
    model,φ,order,cppdegree=cppdegree,trim=trim,limitstol=limitstol)
  node_to_dof_order(cps,fespace,model,order)
end

function compute_closest_point_projections(model::CartesianDiscreteModel,
                                           φ::AlgoimCallLevelSetFunction,
                                           order::Int;
                                           cppdegree::Int=2,
                                           trim::Bool=false,
                                           limitstol::Float64=1.0e-8)
  cdesc = get_cartesian_descriptor(model)
  xmin = cdesc.origin
  xmax = xmin + Point(cdesc.sizes .* cdesc.partition)
  partition = Int32[cdesc.partition...] .* Int32(order)
  fill_cpp_data(φ,partition,xmin,xmax,cppdegree,trim,limitstol,order=order)
end

function compute_closest_point_projections(model::DistributedCartesianDiscreteModel,
                                           φ::AlgoimCallLevelSetFunction{<:Function,<:Function},
                                           order::Int;
                                           cppdegree::Int=2,
                                           trim::Bool=false,
                                           limitstol::Float64=1.0e-8)
  gdesc = model.metadata.descriptor
  xmin = gdesc.origin
  xmax = xmin + Point(gdesc.sizes .* gdesc.partition)
  gpartition = Int32[gdesc.partition...] .* Int32(order)
  map(local_views(model)) do m
    cdesc = get_cartesian_descriptor(m)
    cmin = ( cdesc.origin - gdesc.origin )
    cmin = Int.(round.(cmin.data ./ gdesc.sizes) .* order)
    cpartition = cdesc.partition .* order
    cmax = cmin .+ cpartition
    fill_cpp_data(φ,gpartition,xmin,xmax,cppdegree,trim,limitstol,
              order=order,rmin=Int32[cmin...],rmax=Int32[cmax...])
  end
end

function compute_closest_point_projections(model::DistributedCartesianDiscreteModel,
                                           φ::AlgoimCallLevelSetFunction{<:CellField,<:CellField},
                                           order::Int;
                                           cppdegree::Int=2,
                                           trim::Bool=false,
                                           limitstol::Float64=1.0e-8)
  @notimplemented "The Algoim LS function must be defined from a DistributedCellField"
end

function compute_closest_point_projections(model::DistributedCartesianDiscreteModel,
                                           φ::DistributedAlgoimCallLevelSetFunction,
                                           order::Int;
                                           cppdegree::Int=2,
                                           trim::Bool=false,
                                           limitstol::Float64=1.0e-8)
  gdesc = model.metadata.descriptor
  xmin = gdesc.origin
  xmax = xmin + Point(gdesc.sizes .* gdesc.partition)
  gpartition = gdesc.partition .* order
  gsizes = gdesc.sizes ./ order
  ggrid = vec(map(i->xmin+Point((Tuple(i).-1).*gsizes),CartesianIndices(gpartition.+1)))
  gvals = evaluate(φ.values,ggrid)
  gpartition = Int32[gpartition...]
  map(local_views(model)) do m
    cdesc = get_cartesian_descriptor(m)
    cmin = ( cdesc.origin - gdesc.origin )
    cmin = Int.(round.(cmin.data ./ gdesc.sizes) .* order)
    cpartition = cdesc.partition .* order
    cmax = cmin .+ cpartition
    fill_cpp_data(gvals,gpartition,xmin,xmax,cppdegree,trim,limitstol,
                  rmin=Int32[cmin...],rmax=Int32[cmax...])
  end
end

function node_to_dof_order(ncps,
                           fespace::FESpace,
                           model::CartesianDiscreteModel,
                           order::Int)

  msg = "Is the FE space order the same as the input order?"
  @assert length(ncps) == num_free_dofs(fespace) msg

  D = num_dims(model)
  ncells = num_cells(model)
  cdesc = get_cartesian_descriptor(model)
  partition = cdesc.partition
  rpartition = partition .* order

  orders = tfill(order,Val{D}())
  ones = tfill(1,Val{D}())
  range = CartesianIndices(orders.+1) .- CartesianIndex(ones) # 0-based
  ldof_to_lnode = get_ldof_to_lnode(orders,D)
  node_partition = rpartition .+ 1

  cell_node_ids = lazy_map(1:ncells) do cellid
    cell_ijk = CartesianIndices(partition)[cellid]
    anchor_node = ( cell_ijk .- CartesianIndex(ones) ) .* order .+ CartesianIndex(ones)
    node_ijk = CartesianIndices(rpartition)[anchor_node] 
    range_cis = node_ijk .+ range
    node_ids = LinearIndices(node_partition)[range_cis]
    node_ids[ldof_to_lnode]
  end
  
  cell_dofs_ids = fespace.cell_dofs_ids
  c1 = array_cache(fespace.cell_dofs_ids)
  c2 = array_cache(cell_node_ids)
  dcps = similar(ncps)
  for cellid in 1:ncells
    dof_ids = getindex!(c1,cell_dofs_ids,cellid)
    node_ids = getindex!(c2,cell_node_ids,cellid)
    dcps[dof_ids] = ncps[node_ids]
  end
  dcps
end

function node_to_dof_order(ncps,
                           fespace::DistributedFESpace,
                           model::DistributedCartesianDiscreteModel,
                           order::Int)
  dcps = map(ncps,local_views(fespace),local_views(model)) do cp,fs,m
    node_to_dof_order(cp,fs,m,order)
  end
  PVector(dcps,partition(fespace.gids)) |> consistent! |> wait
  dcps
end

function get_ldof_to_lnode(orders,D)

  # Generate indices of n-faces and order s.t.
  # (1) dimension-increasing (2) lexicographic
  bin_rang_nfaces = tfill(0:1,Val{D}())
  bin_ids_nfaces = vec(collect(Iterators.product(bin_rang_nfaces...)))
  sum_bin_ids_nfaces = sum.(bin_ids_nfaces)
  bin_ids_nfaces = permute!(bin_ids_nfaces,sortperm(sum_bin_ids_nfaces))

  # Generate LIs of basis funs s.t. order by n-faces
  lids_b = LinearIndices(Tuple([orders[i]+1 for i=1:D]))

  eet = eltype(eltype(bin_ids_nfaces))
  f(x) = Tuple( x[i] == one(eet) ? (0:0) : (1:orders[i]:(orders[i]+1)) for i in 1:length(x) )
  g(x) = Tuple( x[i] == one(eet) ? (2:orders[i]) : (0:0) for i in 1:length(x) )
  rang_nfaces = map(f,bin_ids_nfaces)
  rang_own_dofs = map(g,bin_ids_nfaces)

  perm = Int64[]
  for i = 1:length(bin_ids_nfaces)
    cis_nfaces = CartesianIndices(rang_nfaces[i])
    cis_own_dofs = CartesianIndices(rang_own_dofs[i])
    for ci in cis_nfaces
      ci = ci .+ cis_own_dofs
      perm = vcat(perm,reshape(lids_b[ci],length(ci)))
    end
  end

  perm
end

function compute_normal_displacement(
    cps::AbstractVector{<:Point},
    phi::AlgoimCallLevelSetFunction,
    fun,
    dt::Float64,
    Ω::Triangulation)
  disps, _ = compute_normal_displacement!(nothing,cps,phi,fun,dt,Ω)
  disps
end

function compute_normal_displacement(
    cps::AbstractArray,
    phi::Union{AlgoimCallLevelSetFunction,DistributedAlgoimCallLevelSetFunction},
    fun::DistributedCellField,
    dt::Float64,
    Ω::DistributedTriangulation)

  # 1. CPP -> Global cell ID
  model = get_background_model(Ω)
  d_desc = model.metadata.descriptor
  xmin = d_desc.origin
  cell_partition = d_desc.partition
  _cell_ci(x) = CartesianIndex(
    Int.(floor.((x.data.-xmin.data)./d_desc.sizes)).+1)
  _cell_id(x) = LinearIndices(cell_partition)[x]
  cp_gids = map(cps) do cp
    map(_cell_id,map(_cell_ci,cp))
  end
  
  # 2. CPP -> Owner rank ID
  ranks = model.metadata.ranks
  mesh_partition = model.metadata.mesh_partition
  ghost = map(i->true,mesh_partition)
  upartition = uniform_partition(ranks,
                                 mesh_partition,
                                 cell_partition,
                                 ghost)
  cp_owner_ids = find_owner(upartition,cp_gids)

  # 3. Exchange graph
  snd_ids = map(i->sort(unique(i)),cp_owner_ids)
  graph = ExchangeGraph(snd_ids)

  # 4. Exchange local indices
  snd_indices = map(snd_ids,cp_owner_ids) do snd_ids,cp_owner_ids
    snd_indices = map(snd_ids) do s
      findall(i->i==s,cp_owner_ids)
    end
    JaggedArray(snd_indices)
  end

  # 5. Exchange CPPs
  snd_cpps = map(cps,snd_indices) do cps,snd_indices
    map(snd_indices) do indices
      map(Reindex(cps),indices)
    end
  end
  rcv_cpps = allocate_exchange(snd_cpps,graph)
  exchange!(rcv_cpps,snd_cpps,graph) |> wait
  
  # 6. Exchange Global cell IDs
  snd_gids = map(cp_gids,snd_indices) do cp_gids,snd_indices
    map(snd_indices) do indices
      map(Reindex(cp_gids),indices)
    end
  end
  rcv_gids = allocate_exchange(snd_gids,graph)
  exchange!(rcv_gids,snd_gids,graph) |> wait

  # 7. Compute local cell IDs
  gids = get_cell_gids(model)
  rcv_lids = map(global_to_local(gids),rcv_gids) do g_to_l,rcv_gids
    map(Reindex(g_to_l),rcv_gids.data)
  end

  # 8. Evaluate displacements
  normal_phi = normal(phi,Ω)
  snd_disps = map(rcv_lids,
      rcv_cpps,
      local_views(Ω),
      local_views(fun),
      local_views(normal_phi)) do point_to_cell,cpps,Ω,fun,normal_phi
    T = eltype(eltype(eltype(cpps)))
    length(cpps) == 0 && return JaggedArray(zeros(T,0),cpps.ptrs)
    cell_to_points, _ = make_inverse_table(point_to_cell,num_cells(Ω))
    cell_to_xs = lazy_map(Broadcasting(Reindex(cpps.data)),cell_to_points)
    cell_point_xs = CellPoint(cell_to_xs,Ω,PhysicalDomain())
    fun_xs = evaluate(fun,cell_point_xs)
    nΓ_xs = evaluate(normal_phi,cell_point_xs)
    cell_point_disp = lazy_map(Broadcasting(⋅),fun_xs,nΓ_xs)
    cache_vals = array_cache(cell_point_disp)
    cache_ctop = array_cache(cell_to_points)
    disps = zeros(Float64,length(cpps.data))
    for cell in 1:length(cell_to_points)
      pts = getindex!(cache_ctop,cell_to_points,cell)
      vals = getindex!(cache_vals,cell_point_disp,cell)
      for (i,pt) in enumerate(pts)
        val = vals[i]
        disps[pt] = dt * val
      end
    end
    JaggedArray(disps,cpps.ptrs)
  end

  # 9. Exchange displacements
  rgraph = reverse(graph)
  rcv_disps = allocate_exchange(snd_disps,rgraph)
  exchange!(rcv_disps,snd_disps,rgraph) |> wait

  # 10. Merge displacements
  disps = map(rcv_disps,snd_indices) do rcv_disps,snd_indices
    disps = zeros(eltype(eltype(rcv_disps)),length(snd_indices.data))
    map(rcv_disps,snd_indices) do d,s
      disps[s] .= d
    end
    disps
  end

end

function compute_normal_displacement!(
    cache,
    cps::AbstractVector{<:Point},
    phi::AlgoimCallLevelSetFunction,
    fun,
    dt::Float64,
    Ω::Triangulation)
  # Note that cps must be (scalar) DoF-numbered, not lexicographic-numbered
  if isnothing(cache)
    searchmethod = KDTreeSearch()
    cache = _point_to_cell_cache(searchmethod,Ω)
  end
  x_to_cell(x) = _point_to_cell!(cache, x)
  point_to_cell = lazy_map(x_to_cell, cps)
  cell_to_points, _ = make_inverse_table(point_to_cell, num_cells(Ω))
  cell_to_xs = lazy_map(Broadcasting(Reindex(cps)), cell_to_points)
  cell_point_xs = CellPoint(cell_to_xs, Ω, PhysicalDomain())
  fun_xs = evaluate(fun,cell_point_xs)
  nΓ_xs = evaluate(normal(phi,Ω),cell_point_xs)
  cell_point_disp = lazy_map(Broadcasting(⋅),fun_xs,nΓ_xs)
  cache_vals = array_cache(cell_point_disp)
  cache_ctop = array_cache(cell_to_points)
  disps = zeros(Float64,length(cps))
  for cell in 1:length(cell_to_points)
    pts = getindex!(cache_ctop,cell_to_points,cell)
    vals = getindex!(cache_vals,cell_point_disp,cell)
    for (i,pt) in enumerate(pts)
      val = vals[i]
      disps[pt] = dt * val
    end
  end
  disps, cache
end

function compute_normal_displacement!(
    cache,
    uₕ::CellField,
    cpₕ::FEFunction,
    phi::AlgoimCallLevelSetFunction,
    dt::Float64,
    Ω::Triangulation)
  # Note that cps must be (scalar) DoF-numbered, not lexicographic-numbered
  if isnothing(cache)
    searchmethod = KDTreeSearch()
    cache = _point_to_cell_cache(searchmethod,Ω)
  end
  x_to_cell(x) = _point_to_cell!(cache, x)
  cps = map(Point,eachcol(reshape(get_free_dof_values(cpₕ),num_dims(Ω),:)))
  point_to_cell = lazy_map(x_to_cell, cps)
  cell_to_points, _ = make_inverse_table(point_to_cell, num_cells(Ω))
  cell_to_xs = lazy_map(Broadcasting(Reindex(cps)), cell_to_points)
  cell_point_xs = CellPoint(cell_to_xs, Ω, PhysicalDomain())
  fun_xs = evaluate(uₕ,cell_point_xs)
  nΓ_xs = evaluate(normal(phi,Ω),cell_point_xs)
  cell_point_disp = lazy_map(Broadcasting(⋅),fun_xs,nΓ_xs)
  cache_vals = array_cache(cell_point_disp)
  cache_ctop = array_cache(cell_to_points)
  disps = zeros(Float64,length(cps))
  for cell in 1:length(cell_to_points)
    pts = getindex!(cache_ctop,cell_to_points,cell)
    vals = getindex!(cache_vals,cell_point_disp,cell)
    for (i,pt) in enumerate(pts)
      val = vals[i]
      disps[pt] = dt * val
    end
  end
  disps, cache
end

@inline signed_distance(φ::AlgoimCallLevelSetFunction,x,y) = signed_distance(φ.φ,x,y)

@inline signed_distance(φ::Function,x,y) = sign(φ(y))*norm(x-y)

@inline signed_distance(φ::T,x,y) where {T<:Number} = sign(φ)*norm(x-y)

function _compute_signed_distance(
    φ::AlgoimCallLevelSetFunction{<:Function,<:Function},
    cps::Vector{<:Point{N,T}},cos::Array{<:Point{N,T},N}) where {N,T}
  _dist(x,y) = signed_distance(φ.φ,x,y)
  map(_dist,cps,cos)
end

function _compute_signed_distance(
    φ::AlgoimCallLevelSetFunction{<:CellField,<:CellField},
    cps::Vector{<:Point{N,T}},cos::Array{<:Point{N,T},N}) where {N,T}
  φs = get_free_dof_values(φ.φ)
  map(signed_distance,φs,cps,cos)
end

function _compute_signed_distance(
    φ::AlgoimCallLevelSetFunction{<:Function,<:Function},
    cps::Vector{T},cos::Vector{T},ndists::Int,D::Int ) where {T}
  map(1:ndists) do i
    c = view(cos,D*(i-1)+1:D*i)
    p = view(cps,D*(i-1)+1:D*i)
    signed_distance(φ,p,c)
  end
end

function _compute_signed_distance(
    φ::AlgoimCallLevelSetFunction{<:CellField,<:CellField},
    cps::Vector{T},cos::Vector{T},ndists::Int,D::Int ) where {T}
  φs = get_free_dof_values(φ.φ)
  msg = """\n
    Check that FE space for LS function 
    and closest point projections match"""
  @assert length(φs) == ndists msg
  map(1:ndists) do i
    c = view(cos,D*(i-1)+1:D*i)
    p = view(cps,D*(i-1)+1:D*i)
    signed_distance(φs[i],p,c)
  end
end

function compute_distance_fe_function(
    bgmodel::CartesianDiscreteModel,
    fespace::FESpace,
    φ::AlgoimCallLevelSetFunction,
    order::Int;
    cppdegree::Int=2)
  cps = compute_closest_point_projections(
    fespace,φ,order,cppdegree=cppdegree)
  rmodel = refine(bgmodel,order)
  cos = get_node_coordinates(rmodel)
  cos = node_to_dof_order(cos,fespace,bgmodel,order)
  dists = _compute_signed_distance(φ,cps,cos) 
  FEFunction(fespace,dists)
end

function compute_distance_fe_function(
    bgmodel::DistributedCartesianDiscreteModel,
    fespace::FESpace,
    φ::AlgoimCallLevelSetFunction,
    order::Int;
    cppdegree::Int=2)
  cps = compute_closest_point_projections(
    fespace,φ,order,cppdegree=cppdegree)
  _dists = map(cps,local_views(fespace),local_views(bgmodel)) do cp,fs,bg
    rm = refine(bg,order)
    cos = get_node_coordinates(rm)
    cos = node_to_dof_order(cos,fs,bg,order)
    _compute_signed_distance(φ,cp,cos)
  end
  dists = PVector(_dists,partition(fespace.gids)) 
  FEFunction(fespace,dists)
end

function compute_distance_fe_function(
    bgmodel::DistributedCartesianDiscreteModel,
    fespace::FESpace,
    φ::DistributedAlgoimCallLevelSetFunction,
    order::Int;
    cppdegree::Int=2)
  cps = compute_closest_point_projections(
    fespace,φ,order,cppdegree=cppdegree)
  _dists = map(cps,
               local_views(fespace),
               local_views(bgmodel),
               local_views(φ)) do cp,fs,bg,φl
    rm = refine(bg,order)
    cos = get_node_coordinates(rm)
    cos = node_to_dof_order(cos,fs,bg,order)
    _compute_signed_distance(φl,cp,cos)
  end
  dists = PVector(_dists,partition(fespace.gids))
  FEFunction(fespace,dists)
end

function compute_distance_fe_function(
    fespace_scalar_type::FESpace,
    fespace_vector_type::FESpace,
    closest_point_projections::FEFunction,
    φ::AlgoimCallLevelSetFunction)
  model = get_background_model(get_triangulation(fespace_scalar_type))
  D = num_dims(model)
  cos = get_free_dof_values(interpolate(identity,fespace_vector_type))
  cps = get_free_dof_values(closest_point_projections)
  ndists = num_free_dofs(fespace_scalar_type)
  dists = _compute_signed_distance(φ,cps,cos,ndists,D)
  FEFunction(fespace_scalar_type,dists)
end

function apply_mask(υₕ::DistributedSingleFieldFEFunction,mask::Function)
  trian = get_triangulation(υₕ)
  gids  = get_cell_gids(get_background_model(trian))
  map(local_views(υₕ),local_views(trian),local_views(gids)) do lυₕ,t,g
    own_to_local = findall(!iszero,local_to_own(g)[t.tface_to_mface])
    cv = lazy_map(Reindex(get_cell_dof_values(lυₕ)),own_to_local)
    mask(cv)
  end
end

function apply_mask(υₕ::SingleFieldFEFunction,mask::Function)
  cv = get_cell_dof_values(υₕ)
  mask(cv)
end

@inline active_mask(cv) =
  lazy_map(cv) do ccv
    all(ccv .> 0.0) ? false : true
  end

@inline narrow_band_mask(cv) =
  lazy_map(cv) do ccv
    maximum(ccv) * minimum(ccv) < 0.0 ? true : false
  end

function active_triangulation(Ω::Triangulation,
                              φ::SingleFieldFEFunction,
                              Vbg::SingleFieldFESpace,
                              is_a::AbstractArray,
                              δ::Float64)
  φʳ = interpolate_everywhere(φ-δ,Vbg)
  is_aʳ = apply_mask(φʳ,active_mask)
  is_nᵃ = lazy_map((a,aʳ)->a|aʳ,is_a,is_aʳ)
  Triangulation(Ω,is_nᵃ),is_nᵃ
end

function active_triangulation(Ω::DistributedTriangulation,
                              φ::DistributedSingleFieldFEFunction,
                              Vbg::DistributedFESpace,
                              is_a::AbstractArray,
                              δ::Float64)
  φʳ = interpolate_everywhere(φ-δ,Vbg)
  is_aʳ = apply_mask(φʳ,active_mask)
  is_nᵃ = map(local_views(is_a),local_views(is_aʳ)) do is_a,is_aʳ
    lazy_map((a,aʳ)->a|aʳ,is_a,is_aʳ)
  end
  Ωˡ = map((lt,ca)->Triangulation(lt,ca),local_views(Ω),is_nᵃ)
  Mbg = get_background_model(Ω)
  DistributedTriangulation(Ωˡ,Mbg),is_nᵃ
end

function active_triangulation(Ω::DistributedTriangulation,
                              φ::DistributedSingleFieldFEFunction,
                              ηʳ::DistributedSingleFieldFEFunction,
                              Vbg::DistributedFESpace,
                              is_a::AbstractArray,
                              δ::Float64)
  φʳ = interpolate_everywhere(φ-δ,Vbg)
  is_aʳ = apply_mask(φʳ,active_mask)
  is_rʳ = apply_mask(ηʳ,active_mask)
  is_aʳ = map(local_views(is_aʳ),local_views(is_rʳ)) do is_aʳ,is_rʳ
    lazy_map((aʳ,rʳ)->(aʳ&(~rʳ)),is_aʳ,is_rʳ)
  end
  is_nᵃ = map(local_views(is_a),local_views(is_aʳ)) do is_a,is_aʳ
    lazy_map((a,aʳ)->(a|aʳ),is_a,is_aʳ)
  end
  Ωˡ = map((lt,ca)->Triangulation(lt,ca),local_views(Ω),is_nᵃ)
  Mbg = get_background_model(Ω)
  DistributedTriangulation(Ωˡ,Mbg),is_nᵃ
end

function narrow_band_triangulation(Ω::Triangulation,
                                   φ::SingleFieldFEFunction,
                                   Vbg::SingleFieldFESpace,
                                   is_c::AbstractArray,
                                   δ::Float64)
  φʳ = interpolate_everywhere(φ-δ,Vbg)
  φˡ = interpolate_everywhere(φ+δ,Vbg)
  is_cʳ = apply_mask(φʳ,narrow_band_mask)
  is_cˡ = apply_mask(φˡ,narrow_band_mask)    
  is_nᶜ = lazy_map((c,cʳ,cˡ)->c|cʳ|cˡ,is_c,is_cʳ,is_cˡ)
  Triangulation(Ω,is_nᶜ),is_nᶜ
end

function narrow_band_triangulation(Ω::DistributedTriangulation,
                                   φ::DistributedSingleFieldFEFunction,
                                   Vbg::DistributedFESpace,
                                   is_c::AbstractArray,
                                   δ::Float64)
  φʳ = interpolate_everywhere(φ-δ,Vbg)
  φˡ = interpolate_everywhere(φ+δ,Vbg)
  is_cʳ = apply_mask(φʳ,narrow_band_mask)
  is_cˡ = apply_mask(φˡ,narrow_band_mask)    
  is_nᶜ = map(local_views(is_c),
              local_views(is_cʳ),
              local_views(is_cˡ)) do is_c,is_cʳ,is_cˡ
    lazy_map((c,cʳ,cˡ)->c|cʳ|cˡ,is_c,is_cʳ,is_cˡ)
  end
  Ωᶜ = map((lt,ca)->Triangulation(lt,ca),local_views(Ω),is_nᶜ)
  Mbg = get_background_model(Ω)
  DistributedTriangulation(Ωᶜ,Mbg),is_nᶜ
end

function narrow_band_triangulation(Ωⁿ::Triangulation,
    φ::AlgoimCallLevelSetFunction{<:CellField,<:CellField},
    dΓ::Measure,δ::Float64,order::Int)
  
  is_c = is_cell_active(dΓ)
  Ωᶜ = Triangulation(Ωⁿ,is_c)
  
  ncell_to_ls_values = get_cell_dof_values(φ.φ)
  msg = """\n
    Check that the triangulation of the LS function 
    and the old narrow-band triangulation match."""
  @assert length(ncell_to_ls_values) == num_cells(Ωⁿ) msg
  ccell_to_ls_values = lazy_map(Reindex(ncell_to_ls_values),findall(is_c))
  
  ccell_to_vertex_ids = get_cell_node_ids(Ωᶜ)
  bgmodel = get_background_model(Ωᶜ)
  polytope = get_polytopes(bgmodel)
  @assert length(polytope) == 1
  lag_reffe = LagrangianRefFE(Float64,polytope[1],order)
  vertex_dofs = reduce(vcat,get_face_own_dofs(lag_reffe,0))
  cvals = Vector{Float64}(undef,length(vertex_dofs))
  ccell_vertex_to_is_in = map(ccell_to_ls_values) do lsv
                            cvals = abs.(view(lsv,vertex_dofs)) 
                            cvals .< δ
                          end
  
  bgmodel_grid_topology = get_grid_topology(bgmodel)
  bgmodel_vertex_to_cell = get_faces(bgmodel_grid_topology,0,num_dims(bgmodel))
  ccell_vertex_to_cell_around = 
    lazy_map(Broadcasting(Reindex(bgmodel_vertex_to_cell)),ccell_to_vertex_ids)
  is_n = falses(num_cells(bgmodel))
  ccell_to_parent_cell = get_cell_to_parent_cell(get_active_model(Ωᶜ))
  is_n[ccell_to_parent_cell] .= true
  function flag_narrow_band_cells!(v_to_ca,v_to_i) 
    v_to_i && ( is_n[v_to_ca] .= true )
    nothing
  end
  map(Broadcasting(flag_narrow_band_cells!),
        ccell_vertex_to_cell_around,ccell_vertex_to_is_in)

  Triangulation(bgmodel,is_n)
end

abstract type QhullType end

struct DelaunayTrian <: QhullType end
const delaunaytrian = DelaunayTrian()

struct ConvexHull <: QhullType end
const convexhull = ConvexHull()

get_flags(::QhullType) = @abstractmethod
get_flags(::DelaunayTrian) = "qhull d Qt Qbb Qc Qz"
get_flags(::ConvexHull) = "qhull Qt Qc"

get_dimension(::QhullType,dim) = @abstractmethod
get_dimension(::DelaunayTrian,dim) = dim
get_dimension(::ConvexHull,dim) = dim-1

using Gridap.Visualization
import Gridap.Visualization: visualization_data
import Gridap.Visualization: writevtk

function visualization_data(meas::Measure,filename;cellfields=Dict(),qhulltype=DelaunayTrian())
  node_coordinates = collect(Iterators.flatten(meas.quad.cell_point.values))
  grid = _to_grid(node_coordinates,qhulltype)
  ndata = Dict()
  for (k,v) in cellfields
    pts = get_cell_points(meas)
    eval = evaluate(v,pts)
    ndata[k] = collect(Iterators.flatten(eval))
  end
  visualization_data(grid,filename,nodaldata=ndata)
end

function visualization_data(meas::Vector{<:Measure},filename;cellfields=Dict(),qhulltype=DelaunayTrian())
  node_coordinates = vcat(map(m->collect(Iterators.flatten(m.quad.cell_point.values)),meas)...)
  grid = _to_grid(node_coordinates,qhulltype)
  ndata = Dict()
  for (k,v) in cellfields
    pts = map(m->get_cell_points(m),meas)
    eval = map(p->evaluate(v,p),pts)
    ndata[k] = vcat(map(e->collect(Iterators.flatten(e)),eval)...)
  end
  visualization_data(grid,filename,nodaldata=ndata)
end

function _to_grid(node_coordinates::Vector{<:Point{Dp,Tp}},qhulltype) where {Dp,Tp}
  d = get_dimension(qhulltype,Dp)
  connectivity = length(node_coordinates) == 0 ? Matrix{Int32}(undef,0,0) :
    delaunay(reinterpret(node_coordinates),get_flags(qhulltype))[1:(d+1),:]
  cell_node_ids = Table(collect(eachcol(connectivity)))
  reffes = [LagrangianRefFE(Float64,Simplex(Val{d}()),1)]
  cell_types = collect(Fill(Int8(1),length(cell_node_ids)))
  UnstructuredGrid(node_coordinates,cell_node_ids,reffes,cell_types)
end

function visualization_data(meas::DistributedMeasure,
                            filename;
                            cellfields=Dict(),
                            qhulltype=DelaunayTrian())
  trian = meas.trian
  cell_gids = get_cell_gids(trian.model)
  vd = map(
    partition(cell_gids),local_views(meas)) do lindices,meas
    node_coordinates = collect(Iterators.flatten(meas.quad.cell_point.values))
    grid = _to_grid(node_coordinates,qhulltype)
    part = part_id(lindices)
    celldata = Dict{Any,Any}()
    # we do not use "part" since it is likely to be used by the user
    if haskey(celldata,"piece")
      @unreachable "piece is a reserved cell data name"
    end
    celldata["piece"] = fill(part,num_cells(grid))
    # ndata = Dict()
    # for (k,v) in cellfields
    #   pts = get_cell_points(meas)
    #   eval = evaluate(v,pts)
    #   ndata[k] = collect(Iterators.flatten(eval))
    # end
    vd = visualization_data(grid,filename,celldata=celldata) # ,nodaldata=ndata)
    @assert length(vd) == 1
    vd[1]
  end
  [DistributedVisualizationData(vd)]
end

function writevtk(
    arg::DistributedMeasure,args...;
    compress=false,append=true,ascii=false,vtkversion=:default,kwargs...)
  parts=get_parts(arg.trian)
  map(visualization_data(arg,args...;kwargs...)) do visdata
    write_vtk_file(
      parts,visdata.grid,visdata.filebase,celldata=visdata.celldata,nodaldata=visdata.nodaldata,
      compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
    )
  end
end

end # module