module ImplicitGeometriesTests

using Gridap
using Gridap.Geometry
using Gridap.ReferenceFEs
using GridapEmbedded: doughnut
using GridapEmbedded: discretize

const R = 1.2
const r = 0.2

geom = doughnut(R,r)

n = 50
partition = (n,n,n)
desc = CartesianDescriptor(geom.pmin,geom.pmax,partition)
grid = CartesianGrid(desc)

x = get_node_coordinates(grid)
geom_x = discretize(geom,x)

end # module
