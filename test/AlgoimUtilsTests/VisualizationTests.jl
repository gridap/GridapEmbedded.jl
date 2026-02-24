module VisualizationTests

  using CxxWrap
  using Gridap
  using GridapEmbedded
  using Gridap.ReferenceFEs

  const IN = -1
  const OUT = 1
  const CUT = 0

  order = 1
  domain = (-1.1,1.1,-1.1,1.1)

  f(x) = (x[1]*x[1]/(0.5*0.5)+x[2]*x[2]/(0.5*0.5)) - 1.0
  function gradf(x::V) where {V}
      V([2.0*x[1]/(0.5*0.5),2.0*x[2]/(0.5*0.5)])
  end

  n = 14
  degree = 3
  cppdegree = -1

  partition = Int32[n,n]
  bgmodel = CartesianDiscreteModel(domain,partition)
  Ω = Triangulation(bgmodel)

  reffeᵠ = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω,reffeᵠ)

  fₕ = interpolate_everywhere(f,V)
  phi = AlgoimCallLevelSetFunction(fₕ,∇(fₕ))

  squad = Quadrature(algoim,phi,degree,phase=CUT)
  _,dΓ = TriangulationAndMeasure(Ω,squad)

  vquad = Quadrature(algoim,phi,degree,phase=IN)
  _,dΩ = TriangulationAndMeasure(Ω,vquad)

  path = mktempdir()
  writevtk(dΓ,joinpath(path,"res_sur"),cellfields=["f"=>fₕ],qhulltype=convexhull)
  writevtk([dΩ,dΓ],joinpath("res_vol"),cellfields=["f"=>fₕ])

end # module