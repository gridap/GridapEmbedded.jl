module EmbeddedDiscretizationsTests

using Gridap
using Test
using Gridap.Geometry
using GridapEmbedded.CSG
using GridapEmbedded.Interfaces
using GridapEmbedded.LevelSetCutters

const R = 0.7
geom = disk(R)
box = get_metadata(geom)
n = 50
partition = (n,n)

model_quad = CartesianDiscreteModel(box.pmin,box.pmax,partition)
model_tet = simplexify(CartesianDiscreteModel(box.pmin,box.pmax,partition))

reffe = ReferenceFE(lagrangian,Float64,1)

for model in (model_quad, model_tet)

  cutgeom = cut(model,geom)

  Ωact_in = Triangulation(cutgeom,ACTIVE)
  Ωact_in = Triangulation(cutgeom,ACTIVE_IN)
  Ωact_out = Triangulation(cutgeom,ACTIVE_OUT)
  Ω = Triangulation(model)
  Ω_in = Triangulation(cutgeom,PHYSICAL)
  Ω_in = Triangulation(cutgeom,PHYSICAL_IN)
  Ω_out = Triangulation(cutgeom,PHYSICAL_OUT)
  Γ = EmbeddedBoundary(cutgeom)
  Λ_in = GhostSkeleton(cutgeom)
  Λ_in = GhostSkeleton(cutgeom,ACTIVE_IN)
  Λ_out = GhostSkeleton(cutgeom,ACTIVE_OUT)

  test_triangulation(Ωact_in)
  test_triangulation(Ωact_out)
  test_triangulation(Ω_in)
  test_triangulation(Γ)
  test_triangulation(Λ_in)
  test_triangulation(Λ_out)

  dΩ_in = Measure(Ω_in,2)
  dΓ = Measure(Γ,2)
  n_Γ = get_normal_vector(Γ)

  vol = sum( ∫(1)*dΩ_in )
  surf = sum( ∫(1)*dΓ )

  @test abs(pi*R^2 - vol) < 1.0e-3
  @test abs(surf - 2*pi*R) < 1.0e-3

  V_in = FESpace(Ωact_in,reffe,conformity=:H1)
  u(x) = x[1] + x[2]
  u_in = interpolate(u,V_in)
  v_in = FEFunction(V_in,rand(num_free_dofs(V_in)))

  # Check divergence theorem
  a = sum( ∫( ∇(v_in)⋅∇(u_in) )*dΩ_in )
  b = sum( ∫( v_in*n_Γ⋅∇(u_in) )*dΓ )
  @test abs(a-b) < 1.0e-9

  d = mktempdir()
  try
  writevtk(Ω,joinpath(d,"trian"),order=2,cellfields=["v_in"=>v_in])
  writevtk(Γ,joinpath(d,"trian_G"),order=2,cellfields=["v_in"=>v_in,"normal"=>n_Γ])
  writevtk(Λ_in,joinpath(d,"trian_Gg_in"),cellfields=["v_in"=>mean(v_in)])
  writevtk(Λ_out,joinpath(d,"trian_Gg_out"))
  writevtk(Ω_in,joinpath(d,"trian_in"),order=2,cellfields=["v_in"=>v_in])
  #writevtk(cutgeom,joinpath(d,"cutgeom"))
  writevtk(Ωact_in,joinpath(d,"Ωact_in"))
  writevtk(Ωact_out,joinpath(d,"Ωact_out"))
  finally
    rm(d,recursive=true)
  end

end

end # module
