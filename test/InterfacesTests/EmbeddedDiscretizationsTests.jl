module EmbeddedDiscretizationsTests

using Gridap
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters

const R = 0.7
geom = disc(R)
n = 50
partition = (n,n)

model = CartesianDiscreteModel(geom.pmin,geom.pmax,partition)

cutdisc = cut(model,geom)

model_in = DiscreteModel(cutdisc)

model_out = DiscreteModel(cutdisc,OUT)

#writevtk(cutdisc,"cutdisc")
#writevtk(model_in,"model_in")
#writevtk(model_out,"model_out")


end # module
