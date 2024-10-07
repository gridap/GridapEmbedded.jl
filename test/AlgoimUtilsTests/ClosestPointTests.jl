module ClosestPointTests

using Test
using CxxWrap
using Gridap
using GridapEmbedded
using PartitionedArrays
using GridapDistributed

using Gridap.ReferenceFEs
using Gridap.Geometry: get_node_coordinates

const IN = -1
const OUT = 1
const CUT = 0

function run_seq_case(degree)

  dbg = "seq"
  @show dbg

  cells = Int32[32,32]
  
  domain = (-0.61,0.61,0.00,1.22)
  bgmodel = CartesianDiscreteModel(domain,cells)
  Ω = Triangulation(bgmodel)
  
  order = 2
  reffeᵠ = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω,reffeᵠ)
  
  f(x) = 4.0 * ( x[1]*x[1] + x[2]*x[2] ) - 1.0

  _fₕ = interpolate_everywhere(f,V)

  phi = AlgoimCallLevelSetFunction(_fₕ,∇(_fₕ))

  fₕ = compute_distance_fe_function(bgmodel,V,phi,order,cppdegree=degree)
  φₕ = AlgoimCallLevelSetFunction(fₕ,∇(fₕ))

  _degree = order < 3 ? 3 : 2*order
  squad = Quadrature(algoim,φₕ,_degree,phase=CUT)
  s_cell_quad,is_c = CellQuadratureAndActiveMask(bgmodel,squad)

  Γ,dΓ = TriangulationAndMeasure(Ω,s_cell_quad,is_c)

  # writevtk(Γ,"cut_seq_$dbg")
  # writevtk(Ω,"dist_seq_$dbg",cellfields=["LS"=>fₕ])

  # @show ∑( ∫( 1.0 )dΓ )

  # writevtk(dΓ,"cut_seq_quad_$dbg")

  cpps = compute_closest_point_projections(V,φₕ,order)
  nΓ = normal(φₕ,Ω)
  dt = 0.1

  _dₕ = compute_normal_displacement(cpps,φₕ,nΓ,dt,Ω)
  dₕ_fv = PVector(_dₕ,partition(V.gids)) 
  consistent!(dₕ_fv) |> wait
  dₕ = FEFunction(V,dₕ_fv)
  id = parts[1]
  writevtk(Ω,"disps_$id",cellfields=["d"=>dₕ])

end

function run_case(distribute,parts,degree)

  dbg = parts[1]
  @show dbg

  ranks = distribute(LinearIndices((prod(parts),)))
  cells = Int32[8,8]
  
  domain = (-0.61,0.61,0.00,1.22)
  bgmodel = CartesianDiscreteModel(ranks,parts,domain,cells)
  Ω = Triangulation(bgmodel)
  
  order = 2
  reffeᵠ = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω,reffeᵠ)
  
  f(x) = 4.0 * ( x[1]*x[1] + x[2]*x[2] ) - 1.0

  _fₕ = interpolate_everywhere(f,V)

  phi = AlgoimCallLevelSetFunction(_fₕ,∇(_fₕ))

  fₕ = compute_distance_fe_function(bgmodel,V,phi,order,cppdegree=degree)
  φₕ = AlgoimCallLevelSetFunction(fₕ,∇(fₕ))

  _degree = order < 3 ? 3 : 2*order
  squad = Quadrature(algoim,φₕ,_degree,phase=CUT)
  s_cell_quad,is_c = CellQuadratureAndActiveMask(bgmodel,squad)

  Ωᶜ,_ = TriangulationAndMeasure(Ω,s_cell_quad,is_c)

  # writevtk(Γ,"cut_$dbg")
  # writevtk(Ω,"dist_$dbg",cellfields=["LS"=>fₕ])

  # @show ∑( ∫( 1.0 )dΓ )

  # writevtk(dΓ,"cut_quad_$dbg",qhulltype=convexhull)

  cpps = compute_closest_point_projections(V,φₕ,order,
                cppdegree=3,trim=true,limitstol=1.0e-2)
  nΓ = normal(φₕ,Ω)
  dt = 0.1

  reffeᶜ = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
  Vᶜ = TestFESpace(Ωᶜ,reffeᶜ)
  fun(x) = cos(atan(x[2],x[1])-pi) * ∇(_fₕ)(x)
  Ψₕ = interpolate_everywhere(fun,Vᶜ)

  Ωⱽ = get_triangulation(V)
  _dₕ = compute_normal_displacement(cpps,φₕ,Ψₕ,dt,Ωⱽ)
  dₕ_fv = PVector(_dₕ,partition(V.gids)) 
  consistent!(dₕ_fv) |> wait
  dₕ = FEFunction(V,dₕ_fv)
  id = parts[1]
  writevtk(Ω,"disps_$id",cellfields=["d"=>dₕ,"p"=>fₕ])

end

end # module