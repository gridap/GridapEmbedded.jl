module CSG

using Gridap.Helpers
using Test
using AbstractTrees

import AbstractTrees: children
import AbstractTrees: printnode

export Node
export UnaryNode
export Leaf
export replace_data
export Geometry
export get_tree
export get_name
export get_metadata
export replace_metadata
export similar_geometry
export compatible_geometries
export test_geometry
export intersection
export get_geometry
export get_geometry_names

include("Nodes.jl")

include("Geometries.jl")

end # module
