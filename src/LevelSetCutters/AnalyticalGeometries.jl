
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

function cylinder(R;x0=zero(Point{3,eltype(R)}),v=VectorValue(1,0,0),name="cylinder")

  d = v/norm(v)

  function cylinderfun(x)
    _cylinder(x,R,x0,d)
  end

  tree = Leaf((cylinderfun,name,nothing))
  AnalyticalGeometry(tree)

end

@inline function _cylinder(x::Point,R,x0,v)
  w = x-x0
  A = w*v
  H2 = w*w
  B = H2-A*A
  B - R^2
end

function plane(;x0=Point(0,0,0),v=VectorValue(1,0,0),name="plane")

  function planefun(x)
    _plane(x,x0,v)
  end

  tree = Leaf((planefun,name,nothing))
  AnalyticalGeometry(tree)
end

@inline function _plane(x::Point,x0,v)
  w = x-x0
  A = w*v
  A
end

function cube(;L=1,x0=Point(0,0,0),name="cube")

  e1 = VectorValue(1,0,0)
  e2 = VectorValue(0,1,0)
  e3 = VectorValue(0,0,1)

  plane1 = plane(x0=x0-0.5*L*e3,v=-e3)
  plane2 = plane(x0=x0+0.5*L*e3,v=+e3)
  plane3 = plane(x0=x0-0.5*L*e2,v=-e2)
  plane4 = plane(x0=x0+0.5*L*e2,v=+e2)
  plane5 = plane(x0=x0-0.5*L*e1,v=-e1)
  plane6 = plane(x0=x0+0.5*L*e1,v=+e1)
  
  geo12 = intersect(plane1,plane2)
  geo34 = intersect(plane3,plane4)
  geo56 = intersect(plane5,plane6)

  intersect(intersect(geo12,geo34),geo56,name=name)

end


