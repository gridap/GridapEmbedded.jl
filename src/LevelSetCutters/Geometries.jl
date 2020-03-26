
struct AnalyticalGeometry{D,T}
  pmin::Point{D,T}
  pmax::Point{D,T}
  intersection::Bool
  ls_to_function::Vector{<:Function}
  ls_to_name::Vector{String}
end

struct DiscreteGeometry{D,T}
  pmin::Point{D,T}
  pmax::Point{D,T}
  intersection::Bool
  ls_to_point_to_value::Vector{Vector{T}}
  ls_to_name::Vector{String}
end

function DiscreteGeometry(geom::DiscreteGeometry,ls_to_point_to_value)
  DiscreteGeometry(geom.pmin,geom.pmax,geom.intersection,ls_to_point_to_value,geom.ls_to_name)
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

  DiscreteGeometry(g.pmin,g.pmax,g.intersection,ls_to_point_to_value,g.ls_to_name)
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
  ls_to_name = ["doughnut"]

  AnalyticalGeometry(pmin,pmax,intersection,ls_to_function,ls_to_name)
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
  ls_to_name = [ "walls", "inlet", "outlet" ]

  AnalyticalGeometry(pmin,pmax,intersection,ls_to_function,ls_to_name)
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

function olympic_rings(R,r)

  z = zero(R)
  x0 = Point(-(r+r+R),z,z)
  g1 = doughnut(R,r,x0)
  x0 = Point(r+r+R,z,z)
  g2 = doughnut(R,r,x0)
  x0 = Point(z,-R+r,z)
  g3 = doughnut(R,r,x0)
  x0 = Point(2*(r+r+R),-R+r,z)
  g4 = doughnut(R,r,x0)
  x0 = Point(3*(r+r+R),z,z)
  g5 = doughnut(R,r,x0)

  f1 = first(g1.ls_to_function)
  f2 = first(g2.ls_to_function)
  f3 = first(g3.ls_to_function)
  f4 = first(g4.ls_to_function)
  f5 = first(g5.ls_to_function)
  ls_to_function = [f1,f2,f3,f4,f5]
  intersection = false
  ls_to_name = ["ring1","ring2","ring3","ring4","ring5"]

  pmin, pmax = _olympic_rings_box(R,r)

  AnalyticalGeometry(pmin,pmax,intersection,ls_to_function,ls_to_name)
end

function _olympic_rings_box(R,r)
  m = 0.1*r
  A = 2*(R+r)+r+m
  B = (R+r)+m
  C = r+m
  pmin = Point(-A,-B-(R-r),-C)
  pmax = Point(A+2*(r+r+R),B,C)
  pmin, pmax
end

