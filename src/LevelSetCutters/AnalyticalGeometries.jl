
struct AnalyticalGeometry <: CSG.Geometry
  tree::Node
end

get_tree(a::AnalyticalGeometry) = a.tree

similar_geometry(a::AnalyticalGeometry,tree::Node) = AnalyticalGeometry(tree)

compatible_geometries(a::AnalyticalGeometry,b::AnalyticalGeometry) = (a,b)

struct BoundingBox{D,T}
  pmin::Point{D,T}
  pmax::Point{D,T}
end

# Factories

function doughnut(R,r;x0=zero(Point{3,typeof(R)}),name="doughnut")

  box = _doughnut_box(R,r,x0)

  function doughnutfun(x)
    _doughnut_fun(x,R,r,x0)
  end

  tree = Leaf((doughnutfun,name,box))

  AnalyticalGeometry(tree)
end

function _doughnut_box(R,r,x0)
  m = 0.1*r
  A = (R+r)+m
  B = r+m
  pmin = Point(-A,-A,-B) + x0
  pmax = Point(A,A,B) + x0
  BoundingBox(pmin, pmax)
end

@inline function _doughnut_fun(x::Point,R,r,x0)
  _x = x - x0
  (R - sqrt(_x[1]^2+_x[2]^2) )^2 + _x[3]^2 - r^2
end

function sphere(R;x0=zero(Point{3,eltype(R)}),name="sphere")

  function spherefun(x)
    _sphere(x,R,x0)
  end

  box = _sphere_box(R,x0)
  tree = Leaf((spherefun,name,box))
  AnalyticalGeometry(tree)
end

function _sphere_box(R,x0)
  e = 1.01
  pmin = x0 - e*R
  pmax = x0 + e*R
  BoundingBox(pmin, pmax)
end

@inline function _sphere(x::Point,R,x0)
  w = x-x0
  A = w*w - R^2
  A
end

function disk(R;x0=zero(Point{2,eltype(R)}),name="disk")

  function diskfun(x)
    _sphere(x,R,x0)
  end

  box = _sphere_box(R,x0)
  tree = Leaf((diskfun,name,box))
  AnalyticalGeometry(tree)
end


