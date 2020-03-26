
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
  bgmodel::DiscreteModel
  bgcell_to_inoutcut::Vector{Int8}
  subcells_in::SubTriangulation{Dp,T}
  subcells_out::SubTriangulation{Dp,T}
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

function writevtk(st::FacetSubTriangulation,filename::String,celldata=[])
  ug = UnstructuredGrid(st)
  degree = 0
  quad = CellQuadrature(ug,degree)
  dS = integrate(1,ug,quad)

  newcelldata = [
    "normal"=>st.facet_to_normal,
    "bgcell"=>st.facet_to_bgcell,
    "dS"=>dS]

  _celldata = vcat(celldata,newcelldata)

  write_vtk_file(ug,filename,celldata=_celldata)
end

function writevtk(cutdisc::EmbeddedDiscretization,filename::String)

  filename_bg = filename * "_background"
  grid = get_grid(cutdisc.bgmodel)
  write_vtk_file(grid,filename_bg,celldata=["inoutcut"=>cutdisc.bgcell_to_inoutcut])

  filename_in = filename * "_subcells_in"
  writevtk(cutdisc.subcells_in,filename_in)

  filename_out = filename * "_subcells_out"
  writevtk(cutdisc.subcells_out,filename_out)

  filename_subfacets = filename * "_subfacets"
  subfacets, facet_to_tag = merge_facet_sub_triangulations(cutdisc.tag_to_subfacets)
  celldata = ["tag" => facet_to_tag]
  ntags = length(cutdisc.tag_to_subfacets)
  for tag in 1:ntags
    name = cutdisc.tag_to_name[tag]
    data = copy(facet_to_tag)
    mask = facet_to_tag .!= tag
    data[mask] .= 0
    push!(celldata,name=>data)
  end
  writevtk(subfacets,filename_subfacets,celldata)

  nothing
end

function split_in_out(st::SubTriangulation)
 st_in = take_in_or_out(st,IN)
 st_out = take_in_or_out(st,OUT)
 st_in, st_out 
end

function take_in_or_out(st::SubTriangulation{D},in_or_out) where D
  ntcells = 0
  for inoutcut in st.cell_to_inoutcut
    if  inoutcut == in_or_out
      ntcells += 1
    end
  end
  tcell_to_inoutcut = fill(Int8(in_or_out),ntcells)
  tcell_to_bgcell = zeros(Int32,ntcells)
  ntlpoints = D + 1
  tcell_to_points_ptrs = fill(Int32(ntlpoints),ntcells+1)
  length_to_ptrs!(tcell_to_points_ptrs)
  ndata = tcell_to_points_ptrs[end]-1
  tcell_to_points_data = zeros(Int,ndata)
  tcell_to_points = Table(tcell_to_points_data,tcell_to_points_ptrs)
  tcell = 0
  for (cell, inoutcut) in enumerate(st.cell_to_inoutcut)
    if  inoutcut == in_or_out
      tcell += 1
      tcell_to_bgcell[tcell] = st.cell_to_bgcell[cell]
      a = tcell_to_points.ptrs[tcell]-1
      b = st.cell_to_points.ptrs[cell]-1
      for i in 1:ntlpoints
        tcell_to_points.data[a+i] = st.cell_to_points.data[b+i]
      end
    end
  end
  SubTriangulation(
    tcell_to_points,
    tcell_to_inoutcut,
    tcell_to_bgcell,
    st.point_to_coords,
    st.point_to_rcoords)
end

function merge_facet_sub_triangulations(tag_to_subfacets)

  fst = empty(first(tag_to_subfacets))
  facet_to_tag = Int8[]

  for (i,fst_i) in enumerate(tag_to_subfacets)
    facet_to_tag_i = fill(Int8(i),length(fst_i.facet_to_normal))
    append!(fst,fst_i)
    append!(facet_to_tag,facet_to_tag_i)
  end

  fst, facet_to_tag
end

function Base.empty(st::FacetSubTriangulation{Dp,T}) where {Dp,T}

  facet_to_points = Table(Int[],Int32[1,])
  facet_to_normal = Point{Dp,T}[]
  facet_to_bgcell = Int32[]
  point_to_coords = Point{Dp,T}[]
  point_to_rcoords = Point{Dp,T}[]

  FacetSubTriangulation(
    facet_to_points,
    facet_to_normal,
    facet_to_bgcell,
    point_to_coords,
    point_to_rcoords)
end

function Base.append!(a::FacetSubTriangulation{D},b::FacetSubTriangulation{D}) where D

  o = length(a.point_to_coords)

  append!(a.facet_to_normal, b.facet_to_normal)
  append!(a.facet_to_bgcell, b.facet_to_bgcell)
  append!(a.point_to_coords, b.point_to_coords)
  append!(a.point_to_rcoords, b.point_to_rcoords)

  nini = length(a.facet_to_points.data)+1
  append!(a.facet_to_points.data, b.facet_to_points.data)
  nend = length(a.facet_to_points.data)
  append_ptrs!(a.facet_to_points.ptrs,b.facet_to_points.ptrs)
  for i in nini:nend
    a.facet_to_points.data[i] += o
  end

  a
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

