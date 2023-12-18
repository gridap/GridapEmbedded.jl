module DualQuadraturesTests

using Test
using Gridap
using Gridap.ReferenceFEs
using GridapEmbedded
using GridapEmbedded.Interfaces

# 3D funs

u¹(x) = ( (x[1]-0.25) * (x[1]-0.25) +
          (x[2]-0.25) * (x[2]-0.25) + 
          (x[3]-0.25) * (x[3]-0.25) ) / ( 0.20 * 0.20 ) - 1.0

u²(x) = ( (x[1]-0.75) * (x[1]-0.75) +
          (x[2]-0.75) * (x[2]-0.75) + 
          (x[3]-0.75) * (x[3]-0.75) ) / ( 0.22 * 0.22 ) - 1.0

u³(x) = ( (x[1]-0.50) * (x[1]-0.50) +
          (x[2]-0.50) * (x[2]-0.50) + 
          (x[3]-0.50) * (x[3]-0.50) ) / ( 0.40 * 0.40 ) - 1.0

# 2D funs

v¹(x) = ( (x[1]-0.25) * (x[1]-0.25) +
          (x[2]-0.25) * (x[2]-0.25) ) / ( 0.20 * 0.20 ) - 1.0

v²(x) = ( (x[1]-0.75) * (x[1]-0.75) +
          (x[2]-0.75) * (x[2]-0.75) ) / ( 0.22 * 0.22 ) - 1.0

v³(x) = ( (x[1]-0.50) * (x[1]-0.50) +
          (x[2]-0.50) * (x[2]-0.50) ) / ( 0.40 * 0.40 ) - 1.0

function run_case_fe_function(domain,cells,order,degree,phase¹,phase²,u¹,u²)

  model = CartesianDiscreteModel(domain,cells)
  Ω = Triangulation(model)

  reffe = ReferenceFE(lagrangian,Float64,order)
  V = TestFESpace(Ω,reffe)
  U = TrialFESpace(V)

  uₕ¹  = interpolate(u¹,U)
  ∇uₕ¹ = ∇(uₕ¹)

  uₕ²  = interpolate(u²,U)
  ∇uₕ² = ∇(uₕ²)

  phi¹ = AlgoimCallLevelSetFunction(uₕ¹,∇uₕ¹)
  phi² = AlgoimCallLevelSetFunction(uₕ²,∇uₕ²)

  quad = Quadrature(algoim,phi¹,phi²,degree,phase1=phase¹,phase2=phase²)
  dΩ = Measure(Ω,quad,data_domain_style=PhysicalDomain())

  s = ∑(∫(1)dΩ)
  if ( phase¹ == IN ) && ( phase² == OUT )
    @test s ≈ 0.03351032164371535
  elseif ( phase¹ == OUT ) && ( phase² == IN )
    @test s ≈ 0.04460223810211001
  elseif ( phase¹ == IN ) && ( phase² == CUT )
    @test_skip s ≈ 0.11301756191551385 # Different precision in GitHub's CI workflow
  elseif ( phase¹ == OUT ) && ( phase² == CUT )
    @test_skip s ≈ 1.8976021283682756 # Different precision in GitHub's CI workflow
  end

end

domain3D = (0.0,1.1,0.0,1.1,0.0,1.1)

run_case_fe_function(domain3D,(24,24,24),2,5,IN,OUT,u¹,u²)
# run_case_fe_function(domain3D,(24,24,24),2,5,OUT,IN,u¹,u²)

run_case_fe_function(domain3D,(24,24,24),2,5,IN,CUT,u¹,u³)
# run_case_fe_function(domain3D,(24,24,24),2,5,OUT,CUT,u¹,u³)

# domain2D = (0.0,1.1,0.0,1.1)
# domain3D = (0.0,1.1,0.0,1.1,0.0,1.1)
# order  = 2
# phaseA = IN
# phaseB = OUT
# u = u¹
# v = u³

# T = pi*16/25 # intersection 3D
# # T = pi*4/5   # intersection 2D
# # T = pi*(0.2^2+0.22^2) # volumes 2D
# # T = (4/3)*pi*(0.2^3+0.22^3) # volumes 3D
# # T = 4*pi*(0.2^2+0.22^2) # non-intersecting surfaces 3D

# for deg in [3,4,5]
#   @info "Case: degree = $deg"
#   for n in [6,12,24,48]
#     A = run_case_fe_function(domain3D,(n,n,n),2,deg,phaseA,CUT,u¹,u³)
#     B = run_case_fe_function(domain3D,(n,n,n),2,deg,phaseB,CUT,u¹,u³)
#     @show A+B-T
#   end
# end

end # module