function _delaunay(points::Vector{Point{D,T}}) where {D,T}
  n = length(points)
  m = zeros(T,D,n)
  for (i,p) in enumerate(points)
    for (j,pj) in enumerate(p)
      m[j,i] = pj
    end
  end
  cells = delaunay(m)
  [ Vector{Int}(cells[:,k]) for k in 1:size(cells,2)]
end
