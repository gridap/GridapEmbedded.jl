
const IN = -1
const OUT = 1
const INTERFACE = 0
const CUT = 0

struct SubTriangulation{Dp,T}
  cell_to_points::Table{Int,Int32}
  cell_to_inoutcut::Vector{Int8}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dp,T}}
end

struct FacetSubTriangulation{Dp,T}
  facet_to_points::Table{Int,Int32}
  facet_to_normal::Vector{Point{Dp,T}}
  facet_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dp,T}}
  point_to_rcoords::Vector{Point{Dp,T}}
end

struct EmbeddedDiscretization{Dp,T} <: GridapType
  background::Union{DiscreteModel,Grid}
  bgcell_to_inoutcut::Vector{Int8}
  subcells::SubTriangulation{Dp,T}
  tag_to_subfacets::Vector{FacetSubTriangulation{Dp,T}}
  tag_to_name::Vector{String}
end

function UnstructuredGrid(st::SubTriangulation{D}) where D
  reffe = LagrangianRefFE(Float64,Simplex(Val{D}()),1)
  cell_types = fill(Int8(1),length(st.cell_to_points))
  UnstructuredGrid(
    st.point_to_coords,
    st.cell_to_points,
    [reffe,],
    cell_types)
end

function UnstructuredGrid(st::FacetSubTriangulation{Dc}) where Dc
  reffe = LagrangianRefFE(Float64,Simplex(Val{Dc-1}()),1)
  cell_types = fill(Int8(1),length(st.facet_to_points))
  UnstructuredGrid(
    st.point_to_coords,
    st.facet_to_points,
    [reffe,],
    cell_types)
end

function writevtk(st::SubTriangulation,filename::String)
  ug = UnstructuredGrid(st)
  degree = 0
  quad = CellQuadrature(ug,degree)
  dV = integrate(1,ug,quad)
  write_vtk_file(ug,filename,celldata=[
    "inoutcut"=>st.cell_to_inoutcut,
    "bgcell"=>st.cell_to_bgcell,
    "dV"=>dV])
end

function writevtk(st::FacetSubTriangulation,filename::String)
  ug = UnstructuredGrid(st)
  degree = 0
  quad = CellQuadrature(ug,degree)
  dS = integrate(1,ug,quad)
  write_vtk_file(ug,filename,celldata=[
    "normal"=>st.facet_to_normal,
    "bgcell"=>st.facet_to_bgcell,
    "dS"=>dS])
end

function Simplex(p::Polytope)
  D = num_cell_dims(p)
  Simplex(Val{D}())
end

function Simplex(::Val{D}) where D
  extrusion = tfill(TET_AXIS,Val{D}())
  ExtrusionPolytope(extrusion)
end

function Simplex(::Val{2})
  TRI
end

function Simplex(::Val{3})
  TET
end

