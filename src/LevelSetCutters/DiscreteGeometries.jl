
"""
    struct DiscreteGeometry{D,T} <: CSG.Geometry
      tree::Node
      point_to_coords::Vector{Point{D,T}}
    end

## Constructors

    DiscreteGeometry(φh::FEFunction,model::DiscreteModel;name::String="")
    DiscreteGeometry(f::Function,model::DiscreteModel;name::String="")

"""
struct DiscreteGeometry{D,T} <: CSG.Geometry
  tree::Node
  point_to_coords::Vector{Point{D,T}}
end

get_tree(a::DiscreteGeometry) = a.tree

similar_geometry(a::DiscreteGeometry,tree::Node) = DiscreteGeometry(tree,a.point_to_coords)

function compatible_geometries(a::DiscreteGeometry,b::DiscreteGeometry)
  @assert a.point_to_coords === b.point_to_coords || a.point_to_coords == b.point_to_coords
  (a,b)
end

function compatible_geometries(a::DiscreteGeometry,b::AnalyticalGeometry)
  a, discretize(b,a.point_to_coords)
end

function compatible_geometries(b::AnalyticalGeometry,a::DiscreteGeometry)
  discretize(b,a.point_to_coords), a
end

_get_subdomain_grid(model::DiscreteModel) = get_grid(model)
_get_subdomain_grid(trian::Triangulation) = _to_grid_portion(get_grid(trian))
_to_grid_portion(grid::Grid) = grid
_to_grid_portion(grid::Geometry.GridView) = GridPortion(grid.parent,grid.cell_to_parent_cell)

function discretize(a::AnalyticalGeometry,bg::Union{DiscreteModel,Triangulation})
  discretize(a,_get_subdomain_grid(bg))
end

function discretize(a::AnalyticalGeometry,grid::Grid)
  discretize(a,get_node_coordinates(grid))
end

function discretize(a::AnalyticalGeometry,point_to_coords::AbstractArray{<:Point})
  discretize(a,collect1d(point_to_coords))
end

function discretize(a::AnalyticalGeometry,point_to_coords::Vector{<:Point})
  tree = get_tree(a)
  j_to_fun, oid_to_j = _find_unique_leaves(tree)
  j_to_ls = [ fun.(point_to_coords) for fun in j_to_fun ]

  function conversion(data)
    f,name,meta = data
    oid = objectid(f)
    j = oid_to_j[oid]
    ls = j_to_ls[j]
    ls, name, meta
  end

  newtree = replace_data(identity,conversion,tree)
  DiscreteGeometry(newtree,point_to_coords)
end

function _find_unique_leaves(tree)
  i_to_fun = map(n->first(n.data),collect(Leaves(tree)))
  i_to_oid = map(objectid,i_to_fun)
  j_to_oid = unique(i_to_oid)
  j_to_i = collect(Int,indexin(j_to_oid,i_to_oid))
  j_to_fun = i_to_fun[j_to_i]
  oid_to_j = Dict{UInt,Int}( [oid=>j for (j,oid) in enumerate(j_to_oid)] )

  j_to_fun, oid_to_j
end

function DiscreteGeometry(f::Function,bg::Union{DiscreteModel,Triangulation};name::String="")
  discretize(AnalyticalGeometry(f;name=name),bg)
end

function DiscreteGeometry(
  point_to_value::AbstractVector,point_to_coords::AbstractVector;name::String="")
  data = (point_to_value,name,nothing)
  tree = Leaf(data)
  DiscreteGeometry(tree,point_to_coords)
end

# NOTE 1: This assumes that the level set φh is 1st order, i.e that there is a 1-to-1 correspondence
# between nodes in the mesh and dofs in φh. 
# Even if we allowed higher order, the cuts are always linear. Not only it would be a waste
# of time to use higher order, but cuts could actually be wrong.
# This might be developped in the future.
# NOTE 2: When φh has complex-valued DOFs, we define the geometry using the real part of the dofs.
function DiscreteGeometry(
  φh::CellField,bg::Union{DiscreteModel,Triangulation};name::String="")
  grid = _get_subdomain_grid(bg)
  point_to_value = real(get_free_dof_values(φh))
  point_to_coords = collect1d(get_node_coordinates(grid))
  @check length(point_to_value) == length(point_to_coords) """
  The number of free DOFs in φh ($(length(point_to_value))) does not match the number of
  nodes in the background mesh ($(length(point_to_coords))).
  Make sure φh is a 1st-order Lagrangian FEFunction defined on the given background.
  """
  DiscreteGeometry(point_to_value,point_to_coords;name)
end
