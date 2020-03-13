
struct AnalyticalGeometry{D,T}
  pmin::Point{D,T}
  pmax::Point{D,T}
  intersection::Bool
  ls_to_function::Vector{<:Function}
end

struct DiscreteGeometry{D,T}
  pmin::Point{D,T}
  pmax::Point{D,T}
  intersection::Bool
  ls_to_point_to_value::Vector{Vector{T}}
end

function DiscreteGeometry(geom::DiscreteGeometry,ls_to_point_to_value)
  DiscreteGeometry(geom.pmin,geom.pmax,geom.intersection,ls_to_point_to_value)
end

function discretize(geom::AnalyticalGeometry,model::DiscreteModel)
  grid = get_grid(model)
  discretize(geom,grid)
end

function discretize(geom::AnalyticalGeometry,grid::Grid)
  x = get_node_coordinates(grid)
  discretize(geom,x)
end

function discretize(g::AnalyticalGeometry,x::AbstractArray)
  T = eltype(eltype(x))
  npoints = length(x)
  ls_to_point_to_value = [ zeros(T,npoints) for i in 1:length(g.ls_to_function) ]
  for (ls,fun) in enumerate(g.ls_to_function)
    point_to_value = ls_to_point_to_value[ls]
    _discretize!(point_to_value,fun,x)
  end

  DiscreteGeometry(g.pmin,g.pmax,g.intersection,ls_to_point_to_value)
end

function _discretize!(point_to_value,fun,x)

  for point in 1:length(x)
    xi = x[point]
    point_to_value[point] = fun(xi)
  end
end

# Factories

function doughnut(R,r,x0=zero(Point{3,typeof(R)}))

  pmin, pmax = _doughnut_box(R,r,x0)
  intersection = true
  function fun(x)
    _doughnut_fun(x,R,r,x0)
  end
  ls_to_function = [ fun,  ]

  AnalyticalGeometry(pmin,pmax,intersection,ls_to_function)
end

@inline function _doughnut_fun(x::Point,R,r,x0)
  _x = x - x0
  (R - sqrt(_x[1]^2+_x[2]^2) )^2 + _x[3]^2 - r^2
end

function _doughnut_box(R,r,x0)
  m = 0.1*r
  A = (R+r)+m
  B = r+m
  pmin = Point(-A,-A,-B) + x0
  pmax = Point(A,A,B) + x0
  (pmin, pmax)
end

function tube(R,L;x0=zero(Point{3,typeof(R)}),v=VectorValue(1,0,0))

  d = v/norm(v)
  pmin, pmax = _tube_box(R,L,x0,d)
  intersection = true
  walls = x -> _tube_walls(x,R,L,x0,d)
  inlet = x -> _plane(x,x0,-d)
  outlet = x -> _plane(x,x0+L*d,d)
  ls_to_function = [ walls, inlet, outlet ]

  AnalyticalGeometry(pmin,pmax,intersection,ls_to_function)
end

@inline function _tube_walls(x::Point,R,L,x0,v)
  w = x-x0
  A = w*v
  H2 = w*w
  B = H2-A*A
  B - R^2
end

function _tube_box(R,L,x0,v)
  pmin = x0 - R
  pmax = x0 + L*v + R
  pmin, pmax
end

@inline function _plane(x::Point,x0,v)
  w = x-x0
  A = w*v
  A
end

