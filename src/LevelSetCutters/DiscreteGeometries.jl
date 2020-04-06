
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
  i_to_fun = map(n->first(n.data),collect(Leaves(tree)))

  i_to_oid = map(objectid,i_to_fun)
  j_to_oid = unique(i_to_oid)
  j_to_i = collect(Int,indexin(j_to_oid,i_to_oid))
  j_to_fun = i_to_fun[j_to_i]
  j_to_ls = [ fun.(point_to_coords) for fun in j_to_fun ]

  oid_to_j = Dict{UInt,Int}( [oid=>j for (j,oid) in enumerate(j_to_oid)] )

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

