module ClosestPointTests

using Test
using CxxWrap
using Gridap
using GridapEmbedded
using PartitionedArrays
using GridapDistributed

using Gridap.Geometry: get_node_coordinates

const IN = -1
const OUT = 1
const CUT = 0

function run_case(distribute,parts,degree)

  ranks = distribute(LinearIndices((prod(parts),)))
  cells = Int32[32,32,32]
  
  domain = (-1.1,1.1,-1.1,1.1,-1.1,1.1)
  bgmodel = CartesianDiscreteModel(ranks,parts,domain,cells)
  Ω = Triangulation(bgmodel)
  
  reffeᵠ = ReferenceFE(lagrangian,Float64,1)
  V = TestFESpace(Ω,reffeᵠ,conformity=:L2)
  
  f(x) = √(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]) - 0.25
  ∇f(x) = ∇(f)(x)
  fₕ = interpolate_everywhere(f,V)
  phi = AlgoimCallLevelSetFunction(f,∇f)
  
  coords = map(local_views(bgmodel)) do t
    _p = reduce(vcat,collect(get_node_coordinates(t)))
    [ Point(_p[i]...) for i in 1:length(_p) ]
  end
  _cpps = map(coords) do c
    [ p-f(p)*∇f(p) for p in c ]
  end
  cpps = compute_closest_point_projections(bgmodel,phi,cppdegree=degree)
  # map(coords,cpps,_cpps,ranks) do c,p,_p,i
  #   writevtk(c,"cpps_$i",nodaldata=["f"=>f.(p),"cp1"=>p,"cp2"=>_p])
  # end

  partials = map(cpps) do c
    maximum(f.(c))
  end
  maxf = reduce(max,partials,init=typemin(Float64))
  partials = map(cpps) do c
    minimum(f.(c))
  end
  minf = reduce(min,partials,init=typemax(Float64))
  @show maxf,minf

end

end # module