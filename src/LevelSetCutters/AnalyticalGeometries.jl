
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

function AnalyticalGeometry(f::Function)
  tree = Leaf((f,string(nameof(f)),nothing))
  AnalyticalGeometry(tree)
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

function popcorn(;
  r0=0.6,
  σ=0.2,
  A=2,
  x0=zero(Point{3,typeof(r0)}),
  name="popcorn")

  box = _popcorn_box(x0,r0)

  function popcornfun(x)
    _popcorn_fun(x,x0,r0,σ,A)
  end

  tree = Leaf((popcornfun,name,box))

  AnalyticalGeometry(tree)
end

function _popcorn_box(x0,R)
  e = 1.5
  pmin = x0 - e*R
  pmax = x0 + e*R
  BoundingBox(pmin, pmax)
end

@inline function _popcorn_fun(_x,x0,r0,σ,A)
  function point_k(k,r0)
    if 0 <= k && k<=4
      α = 2*k*π/5
      (r0/sqrt(5))*Point(2*cos(α),2*sin(α),1.)
    elseif 5<=k && k <=9
      α = (2*(k-5)-1)*π/5
      (r0/sqrt(5))*Point(2*cos(α),2*sin(α),-1.)
    elseif k==10
      Point(0.,0.,r0)
    else
      Point(0.,0.,-r0)
    end
  end
  x,y,z = _x - x0
  val = sqrt(x^2+y^2+z^2) - r0
  for k in 0:11
    xk,yk,zk = point_k(k,r0)
    val -= A*exp(-((x-xk)^2+(y-yk)^2+(z-zk)^2)/σ^2)
  end
  val
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
  A = w⋅w - R^2
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
  A = w⋅v
  H2 = w⋅w
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
  A = w⋅v
  A
end

function square(;L=1,x0=Point(0,0),name="square",edges=["edge$i" for i in 1:4])

  e1 = VectorValue(1,0)
  e2 = VectorValue(0,1)

  plane1=plane(x0=x0-0.5*L*e2,v=-e2,name=edges[1])
  plane2=plane(x0=x0+0.5*L*e2,v=+e2,name=edges[2])
  plane3=plane(x0=x0-0.5*L*e1,v=-e1,name=edges[3])
  plane4=plane(x0=x0+0.5*L*e1,v=+e1,name=edges[4])

  geo12 = intersect(plane1,plane2)
  geo34 = intersect(plane3,plane4)

  intersect(geo12,geo34,name=name)

end

function quadrilateral(;x0=Point(0,0),d1=VectorValue(1,0),d2=VectorValue(0,1),name="quadrilateral")

    x1 = x0+d1
    x2 = x0+d2
    slope1 = d1[2]/d1[1]
    slope2 = d2[2]/d2[1]

    if slope1 > slope2
      temp = slope1
      slope1 = slope2
      slope2 =temp
      var = x1
      x1 = x2
      x2 = var
    end

    slope_n1 = -1/slope1
    slope_n2 = -1/slope2

    den1 = sqrt(1+slope_n1*slope_n1)
    den2 = sqrt(1+slope_n2*slope_n2)

    n1 = VectorValue(1/den1,slope_n1/den1)
    n2 = VectorValue(1/den2,slope_n2/den2)

    if slope_n1 == -Inf
      n1 = VectorValue(0.0,-1.0)
    end

    if slope_n2 == -Inf
      n2 = VectorValue(0.0,-1.0)
    end

    plane1=plane(x0=x0,v=+n1,name="edge1")
    plane2=plane(x0=x2,v=-n1,name="edge2")
    plane3=plane(x0=x0,v=-n2,name="edge3")
    plane4=plane(x0=x1,v=+n2,name="edge4")

    geo12 = intersect(plane1,plane2)
    geo34 = intersect(plane3,plane4)

    intersect(geo12,geo34,name=name)

end

function cube(;L=1,x0=Point(0,0,0),name="cube")

  e1 = VectorValue(1,0,0)
  e2 = VectorValue(0,1,0)
  e3 = VectorValue(0,0,1)

  plane1 = plane(x0=x0-0.5*L*e3,v=-e3,name="face1")
  plane2 = plane(x0=x0+0.5*L*e3,v=+e3,name="face2")
  plane3 = plane(x0=x0-0.5*L*e2,v=-e2,name="face3")
  plane4 = plane(x0=x0+0.5*L*e2,v=+e2,name="face4")
  plane5 = plane(x0=x0-0.5*L*e1,v=-e1,name="face5")
  plane6 = plane(x0=x0+0.5*L*e1,v=+e1,name="face6")

  geo12 = intersect(plane1,plane2)
  geo34 = intersect(plane3,plane4)
  geo56 = intersect(plane5,plane6)

  intersect(intersect(geo12,geo34),geo56,name=name)

end

function tube(R,L;x0=zero(Point{3,typeof(R)}),v=VectorValue(1,0,0),name="tube")

  d = v/norm(v)
  box = _tube_box(R,L,x0,d)

  walls = cylinder(R,x0=x0,v=d,name="walls")
  inlet = plane(x0=x0,v=-d,name="inlet")
  outlet = plane(x0=x0+L*d,v=d,name="outlet")
  intersect(intersect(inlet,outlet),walls,name=name,meta=box)

end

function _tube_box(R,L,x0,v)
  pmin = x0 - R
  pmax = x0 + L*v + R
  BoundingBox(pmin, pmax)
end

function olympic_rings(R,r,name="olympic_rings")

  box = _olympic_rings_box(R,r)

  z = zero(R)
  geo1 = doughnut(R,r,name="ring1",x0=Point(-(r+r+R),z,z))
  geo2 = doughnut(R,r,name="ring2",x0=Point(r+r+R,z,z))
  geo3 = doughnut(R,r,name="ring3",x0=Point(z,-R+r,z))
  geo4 = doughnut(R,r,name="ring4",x0=Point(2*(r+r+R),-R+r,z))
  geo5 = doughnut(R,r,name="ring5",x0=Point(3*(r+r+R),z,z))

  geo12 = union(geo1,geo2)
  geo123 = union(geo12,geo3)
  geo1234 = union(geo123,geo4)
  geo12345 = union(geo1234,geo5,name=name,meta=box)

  geo12345

end

function _olympic_rings_box(R,r)
  m = 0.1*r
  A = 2*(R+r)+r+m
  B = (R+r)+m
  C = r+m
  pmin = Point(-A,-B-(R-r),-C)
  pmax = Point(A+2*(r+r+R),B,C)
  BoundingBox(pmin, pmax)
end
