
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

