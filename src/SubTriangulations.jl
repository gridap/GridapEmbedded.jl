
struct SubCell
  num_points::Int
  edge_to_points::Vector{Vector{Int}}
end

struct SubTriangulation{Dc,T}
  cell_to_points::Table{Int,Int32}
  cell_to_inoutcut::Vector{Int8}
  cell_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dc,T}}
end

struct FacetSubTriangulation{Dc,T}
  facet_to_points::Table{Int,Int32}
  facet_to_normal::Vector{Point{Dc,T}}
  facet_to_bgcell::Vector{Int32}
  point_to_coords::Vector{Point{Dc,T}}
end
