module GeometryTests

using Test
using AbstractTrees
using GridapEmbedded.CSG
import GridapEmbedded.CSG: get_tree, get_metadata, compatible_geometries, similar_geometry

struct MockGeometry <: Geometry
  tree::Node
end

function MockGeometry(name::String)
  tree = Leaf((name,name,time()))
  MockGeometry(tree)
end

get_tree(a::MockGeometry) = a.tree
similar_geometry(a::MockGeometry,tree::Node) = MockGeometry(tree)
compatible_geometries(a::MockGeometry,b::MockGeometry) = (a,b)

a = MockGeometry("a")
b = MockGeometry("b")
c = union(a,b,name="c")
t = time()
_d = "d"
d = intersect(a,c,name=_d,meta=t)
test_geometry(d)
@test _d === get_name(d)
@test t === get_metadata(d)

#print_tree(stdout,get_tree(d))

meta = "hi!"
e = replace_metadata(d,meta)
test_geometry(e)
@test get_metadata(e) === meta

#print_tree(stdout,get_tree(e))

t = time()
f = !(d,name="f",meta=t)
test_geometry(f)
#print_tree(stdout,get_tree(f))

meta = "hi!"
g = replace_metadata(f,meta)
test_geometry(g)
@test get_metadata(g) === meta

end # module
