
abstract type Geometry end

function get_tree(geo::Geometry)
  @abstractmethod
end

function similar_geometry(a::Geometry,tree::Node)
  @abstractmethod
end

function compatible_geometries(a::Geometry,b::Geometry)
  @abstractmethod
end

function test_geometry(a::Geometry)
  @test isa(get_tree(a), Node)
  get_metadata(a)
  geo = similar_geometry(a,get_tree(a))
  @test isa(geo,Geometry)
  _a,_b = compatible_geometries(a,a)
end

function get_name(geo::Geometry)
  node = get_tree(geo)
  _,name, = node.data
  name
end

function get_metadata(geo::Geometry)
  node = get_tree(geo)
  _,_,meta = node.data
  meta
end

function replace_metadata(geo::Geometry,meta)
  tree = get_tree(geo)
  new_tree = _replace_metadata(tree,meta)
  similar_geometry(geo,new_tree)
end

function _replace_metadata(tree::Node,meta)
  g,name, = tree.data
  data = (g,name,meta)
  Node(data,tree.leftchild,tree.rightchild)
end

function _replace_metadata(tree::Leaf,meta)
  g,name, = tree.data
  data = (g,name,meta)
  Leaf(data)
end

function Base.union(a::Geometry,b::Geometry;name::String="",meta=nothing)
  _a, _b = compatible_geometries(a,b)
  tree = union(get_tree(_a),get_tree(_b),name,meta)
  similar_geometry(_a,tree)
end

function Base.intersect(a::Geometry,b::Geometry;name::String="",meta=nothing)
  _a, _b = compatible_geometries(a,b)
  tree = intersect(get_tree(_a),get_tree(_b),name,meta)
  similar_geometry(_a,tree)
end

function Base.setdiff(a::Geometry,b::Geometry;name::String="",meta=nothing)
  _a, _b = compatible_geometries(a,b)
  tree = setdiff(get_tree(_a),get_tree(_b),name,meta)
  similar_geometry(_a,tree)
end

const intersection = intersect

function Base.union(a::Node,b::Node,name::String,meta)
  Node((:∪, name, meta),a,b)
end

function Base.intersect(a::Node,b::Node,name::String,meta)
  Node((:∩, name, meta),a,b)
end

function Base.setdiff(a::Node,b::Node,name::String,meta)
  Node((:-, name, meta),a,b)
end

