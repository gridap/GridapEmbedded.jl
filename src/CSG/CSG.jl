module CSG

using Gridap.Helpers
using Test

import AbstractTrees: children
import AbstractTrees: printnode

export Node
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

include("Nodes.jl")

include("Geometries.jl")

end # module
