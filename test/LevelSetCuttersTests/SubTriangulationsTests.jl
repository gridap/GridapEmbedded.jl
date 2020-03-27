module SubTriangulationsTests

using Gridap
using Gridap.Geometry
using Gridap.Visualization
using GridapEmbedded.LevelSetCutters
using GridapEmbedded.LevelSetCutters: initial_sub_triangulation
using GridapEmbedded.LevelSetCutters: cut_sub_triangulation

#const R = 1.2
#const r = 0.2
#geom = doughnut(R,r)
#n = 25
#partition = (2*n,2*n,n)

const R = 0.7
const L = 5
geom = tube(R,L,x0=Point(-0.5,0.0,-0.25),v=VectorValue(2,1,1))
n = 50
partition = (n,n,n)

#const R = 1.2
#const r = 0.2
#n = 20
#partition =(8*n,4*n,n)
#geom = olympic_rings(R,r)

grid = CartesianGrid(geom.pmin,geom.pmax,partition)

subtrian, subgeom, bgcell_to_inoutcut = initial_sub_triangulation(grid,geom)

st, ls_to_fst = cut_sub_triangulation(subtrian,subgeom)

write_vtk_file(grid,"grid")
writevtk(st,"st")

for (i,fst) in enumerate(ls_to_fst)
  writevtk(fst,"fst_$i")
end

end #module
