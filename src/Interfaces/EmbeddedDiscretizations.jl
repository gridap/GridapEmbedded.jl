
struct EmbeddedDiscretization{Dp,T} <: GridapType
  bgmodel::DiscreteModel
  bgcell_to_inoutcut::Vector{Int8}
  subcells_in::SubTriangulation{Dp,T}
  subcells_out::SubTriangulation{Dp,T}
  tag_to_subfacets::Vector{FacetSubTriangulation{Dp,T}}
  tag_to_name::Vector{String}
end

function DiscreteModel(cut::EmbeddedDiscretization)
  DiscreteModel(cut,IN)
end

function DiscreteModel(cut::EmbeddedDiscretization,in_or_out)
  pred = i-> (i==CUT) || i==in_or_out
  cell_list = findall(pred, cut.bgcell_to_inoutcut)
  DiscreteModel(cut.bgmodel,cell_list)
end

function GhostSkeleton(cut::EmbeddedDiscretization)
  GhostSkeleton(cut,IN)
end

function GhostSkeleton(cut::EmbeddedDiscretization,in_or_out)

  @assert in_or_out in (IN,OUT)
  cell_to_inoutcut = cut.bgcell_to_inoutcut
  model = cut.bgmodel
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  facet_to_cells = Table(get_faces(topo,D-1,D))
  facet_to_mask = fill(false,length(facet_to_cells))
  _fill_ghost_skeleton_mask!(facet_to_mask,facet_to_cells,cell_to_inoutcut,in_or_out)

  SkeletonTriangulation(model,facet_to_mask)
end

function _fill_ghost_skeleton_mask!(facet_to_mask,facet_to_cells::Table,cell_to_inoutcut,in_or_out)

  nfacets = length(facet_to_cells)
  for facet in 1:nfacets
    a = facet_to_cells.ptrs[facet]
    b = facet_to_cells.ptrs[facet+1]
    ncells_around = b-a
    ncells_around_cut = 0
    ncells_around_active = 0
    for cell_around in 1:ncells_around
      cell = facet_to_cells.data[a-1+cell_around]
      inoutcut = cell_to_inoutcut[cell]
      if (inoutcut == CUT)
        ncells_around_cut += 1
      end
      if (inoutcut == CUT) || (inoutcut == in_or_out)
        ncells_around_active += 1
      end
    end
    if (ncells_around_cut >0) && (ncells_around_active == 2)
      facet_to_mask[facet] = true
    end
  end

end

function Triangulation(cut::EmbeddedDiscretization)
  Triangulation(cut,IN)
end

function Triangulation(cut::EmbeddedDiscretization,in_or_out)
  if in_or_out == IN
    st = cut.subcells_in
  else
    @assert in_or_out == OUT
    st = cut.subcells_out
  end
  trian_cut = SubTriangulationWrapper(st)
  trian = Triangulation(cut.bgmodel)
  cell_to_mask = collect(Bool,cut.bgcell_to_inoutcut .== in_or_out)
  trian_in_or_out = RestrictedTriangulation(trian,cell_to_mask)
  lazy_append(trian_cut,trian_in_or_out)
end

function EmbeddedBoundary(cut::EmbeddedDiscretization,name::String)
  tag = findfirst(i->i==name,cut.tag_to_name)
  EmbeddedBoundary(cut,tag)
end

function EmbeddedBoundary(cut::EmbeddedDiscretization,names::Vector{String})
  tags = findall(i->i in names,cut.tag_to_name)
  EmbeddedBoundary(cut,tags)
end

function EmbeddedBoundary(cut::EmbeddedDiscretization,tags::Vector{<:Integer})
  if length(tags) == 1
    tag = first(tags)
    return EmbeddedBoundary(cut,tag)
  else
    fst = empty(first(cut.tag_to_subfacets))
    for (tag_i,fst_i) in enumerate(cut.tag_to_subfacets)
      if tag_i in tags
        append!(fst,fst_i)
      end
    end
    return FacetSubTriangulationWrapper(fst)
  end
end

function EmbeddedBoundary(cut::EmbeddedDiscretization,tag::Integer)
  FacetSubTriangulationWrapper(cut.tag_to_subfacets[tag])
end

function EmbeddedBoundary(cut::EmbeddedDiscretization)
  EmbeddedBoundary(cut,cut.tag_to_name)
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

function cell_measure(trian_Ω1,n_bgcells)
  trian_cut_1 = trian_Ω1
  quad_cut_1 = CellQuadrature(trian_cut_1,0)
  subcell1_to_dV = integrate(1,trian_cut_1,quad_cut_1)
  subcell1_to_bgcell = get_cell_id(trian_cut_1)
  bgcell_to_dV = zeros(n_bgcells)
  _meas_K_fill!(bgcell_to_dV,subcell1_to_dV,subcell1_to_bgcell)
  bgcell_to_dV
end

function _meas_K_fill!(bgcell_to_dV,subcell1_to_dV,subcell1_to_bgcell)
  for (subcell1, bgcell) in enumerate(subcell1_to_bgcell)
    bgcell_to_dV[bgcell] += subcell1_to_dV[subcell1]
  end
end
