
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

function discretize(a::AnalyticalGeometry,model::DiscreteModel)
  discretize(a,get_grid(model))
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

function _get_value_at_coords(φh::CellField,model::DiscreteModel{Dc,Dp}) where {Dc,Dp}
  @assert DomainStyle(φh) == ReferenceDomain()
  # Cell-to-node map for the original model
  c2n_map = collect1d(get_cell_node_ids(model))

  # Cell-wise node coordinates (in ReferenceDomain coordinates)
  cell_reffe = get_cell_reffe(model)
  cell_node_coords = lazy_map(get_node_coordinates,cell_reffe)

  # Get cell data
  φh_data = CellData.get_data(φh)
  T = return_type(testitem(CellData.get_data(φh)),testitem(testitem(cell_node_coords)))
  values  = Vector{T}(undef,num_nodes(model))
  cell_node_coords_cache = array_cache(cell_node_coords)
  # Loop over cells
  for cell in eachindex(c2n_map)
    field = φh_data[cell]
    node_coords = getindex!(cell_node_coords_cache,cell_node_coords,cell)
    for (iN,node) in enumerate(c2n_map[cell])
      values[node] = field(node_coords[iN])
    end
  end
  return values
end

function DiscreteGeometry(
  point_to_value::AbstractVector,point_to_coords::AbstractVector;name::String="")
  data = (point_to_value,name,nothing)
  tree = Leaf(data)
  DiscreteGeometry(tree,point_to_coords)
end

function DiscreteGeometry(
  φh::CellField,model::DiscreteModel;name::String="")
  point_to_value = _get_value_at_coords(φh,model)
  point_to_coords = collect1d(get_node_coordinates(model))
  DiscreteGeometry(point_to_value,point_to_coords;name)
end