module LevelSetCuttersTests

using Gridap
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters

const R = 0.7
const L = 5
geom = tube(R,L,x0=Point(-0.5,0.0,-0.25),v=VectorValue(2,1,1))
n = 50
partition = (n,n,n)

model = CartesianDiscreteModel(geom.pmin,geom.pmax,partition)

cutter = LevelSetCutter()

cutdisc = cut(cutter,model,geom)

writevtk(cutdisc,"cutdisc")

end # module
