module LookupTablesTests

using MiniQhull
import MiniQhull: delaunay
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.Helpers

function num_cases(nvertices)
  2^nvertices
end

function isout(v)
  v > 0
end

function compute_case(values)
  case = 1
  for (i,v) in enumerate(values)
    if  isout(v)
      case += 2^(i-1)
    end
  end
  case
end

struct LookupTable
  subcells::Vector{Vector{Vector{Int}}}
  inoutcells::Vector{Vector{Bool}}
  subfaces::Vector{Vector{Vector{Int}}}
end

function LookupTable(p::Polytope)

  nvertices = num_vertices(p)

  edge_to_vertices = get_faces(p,1,0)
  vertex_to_coords = get_vertex_coordinates(p)

  cis = CartesianIndices(Tuple(fill(2,nvertices)))
  for ci in cis
    vertex_to_value = 2 .* (Tuple(ci) .-1 ) .- 1

    points = _compute_delaunay_points(vertex_to_value, vertex_to_coords, edge_to_vertices)

  end

end

function delaunay(points::Vector{Point{D,T}}) where {D,T}
  n = length(points)
  m = zeros(T,D,n)
  for (i,p) in enumerate(points)
    for (j,pj) in enumerate(p)
      m[j,i] = pj
    end
  end
  cells = delaunay(m)
  [ cells[:,k] for k in 1:size(cells,2)]
end

function _compute_delaunay_points(vertex_to_value, vertex_to_coords, edge_to_vertices)

  point_to_coords = copy(vertex_to_coords)
  for vertices in edge_to_vertices
    v1 = vertex_to_value[vertices[1]]
    v2 = vertex_to_value[vertices[2]]
    if isout(v1) ⊻ isout(v2)
      p1 = vertex_to_coords[vertices[1]]
      p2 = vertex_to_coords[vertices[2]]
      p = 0.5*(p1+p2)
      push!(point_to_coords,p)
    end
  end

  point_to_coords
end



using Test

@test num_cases(3) == 8
@test compute_case([0,0,0]) == 1
@test compute_case([1,0,0]) == 2
@test compute_case([1,0,1]) == 6
@test compute_case([1,1,1]) == 8

LookupTable(TRI)



end # module
