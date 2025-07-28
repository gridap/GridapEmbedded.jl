
"""
    abstract type Geometry

Abstract type for the definition of a geometry.

## Interface

- `get_tree(geo::Geometry)`
- `similar_geometry(a::Geometry,tree::Node)`
- `compatible_geometries(a::Geometry,b::Geometry)`

"""
abstract type Geometry end

"""
    get_tree(geo::Geometry)
"""
function get_tree(geo::Geometry)
  @abstractmethod
end

"""
    similar_geometry(a::Geometry,tree::Node)
"""
function similar_geometry(a::Geometry,tree::Node)
  @abstractmethod
end

"""
    compatible_geometries(a::Geometry,b::Geometry)
"""
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

"""
    Base.union(a::Geometry,b::Geometry;name::String="",meta=nothing)
"""
function Base.union(a::Geometry,b::Geometry;name::String="",meta=nothing)
  _a, _b = compatible_geometries(a,b)
  tree = union(get_tree(_a),get_tree(_b),name,meta)
  similar_geometry(_a,tree)
end

"""
    Base.intersect(a::Geometry,b::Geometry;name::String="",meta=nothing)
"""
function Base.intersect(a::Geometry,b::Geometry;name::String="",meta=nothing)
  _a, _b = compatible_geometries(a,b)
  tree = intersect(get_tree(_a),get_tree(_b),name,meta)
  similar_geometry(_a,tree)
end

"""
    Base.setdiff(a::Geometry,b::Geometry;name::String="",meta=nothing)
"""
function Base.setdiff(a::Geometry,b::Geometry;name::String="",meta=nothing)
  _a, _b = compatible_geometries(a,b)
  tree = setdiff(get_tree(_a),get_tree(_b),name,meta)
  similar_geometry(_a,tree)
end

"""
    Base.:!(a::Geometry;name::String="",meta=nothing)
"""
function Base.:!(a::Geometry;name::String="",meta=nothing)
  tree = !(get_tree(a),name,meta)
  similar_geometry(a,tree)
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

function Base.:!(a::Node,name::String,meta)
  Node((:!, name, meta),a)
end

function get_geometry(a::Geometry,name::String)
  tree = get_tree(a)
  treei = get_geometry_node(tree,name)
  similar_geometry(a,treei)
end

function get_geometry_names(geo::Geometry)
  tree = get_tree(geo)
  names = String[]
  for a in PreOrderDFS(tree)
    _, name, = a.data
    if name != ""
      push!(names,name)
    end
  end
  names
end

function get_geometry_node(a::Node,name::String)
  for ai in PreOrderDFS(a)
    _, namei, = ai.data
    if namei == name
      return ai
    end
  end
  @unreachable "There is no entity called $name"
end
