module ImplicitGeometriesTests

using Gridap
using Gridap.Geometry
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.LevelSetCutters: discretize

const R = 1.2
const r = 0.2

geom = doughnut(R,r)

n = 50
partition = (n,n,n)
grid = CartesianGrid(geom.pmin,geom.pmax,partition)

geom_x = discretize(geom,grid)

end # module
