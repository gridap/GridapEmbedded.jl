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
  cells = Int32[16,16]
  
  domain = (-1.1,1.1,-1.1,1.1)
  bgmodel = CartesianDiscreteModel(ranks,parts,domain,cells)
  Ω = Triangulation(bgmodel)
  
  order = 2
  reffeᵠ = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω,reffeᵠ)
  
  f(x) = x[1]*x[1]+x[2]*x[2] - 0.25
  # f(x) = x[1]*x[1]+x[2]*x[2]+x[3]*x[3] - 0.25
  ∇f(x) = ∇(f)(x)
  phi = AlgoimCallLevelSetFunction(f,∇f)

  fₕ = compute_distance_fe_function(bgmodel,V,phi,order,cppdegree=degree)
  # id = parts[1]
  # writevtk(Ω,"distance",cellfields=["f"=>fₕ])
  φₕ = AlgoimCallLevelSetFunction(fₕ,∇(fₕ))

  gₕ = compute_distance_fe_function(bgmodel,V,φₕ,order,cppdegree=degree)
  # id = parts[1]
  # writevtk(Ω,"distance_$id",cellfields=["g"=>gₕ])

  # cpps = compute_closest_point_projections(bgmodel,φₕ)
  cpps = compute_closest_point_projections(V,φₕ,order)
  # cpps = compute_closest_point_projections(V,phi,order)
  # map(cpps) do cpp
  #   @show cpp
  # end
  nΓ = normal(φₕ,Ω)
  dt = 0.1

  _dₕ = compute_normal_displacement(cpps,φₕ,nΓ,dt,Ω)
  dₕ_fv = PVector(_dₕ,partition(V.gids)) 
  consistent!(dₕ_fv) |> wait
  dₕ = FEFunction(V,dₕ_fv)
  id = parts[1]
  writevtk(Ω,"disps_$id",cellfields=["d"=>dₕ])

  # dₕ = compute_normal_displacement(cpps,phi,nΓ,dt,Ω)
  # map(dₕ) do d
  #   @show d
  # end

end

end # module