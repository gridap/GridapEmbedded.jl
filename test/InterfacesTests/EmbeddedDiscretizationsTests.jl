module EmbeddedDiscretizationsTests

using Gridap
using Test
using Gridap.Arrays
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
  @test isa(Ωact_in,Gridap.Geometry.BodyFittedTriangulation)
  Ωact_in = Triangulation(cutgeom,ACTIVE_IN)
  @test isa(Ωact_in,Gridap.Geometry.BodyFittedTriangulation)
  Ωact_out = Triangulation(cutgeom,ACTIVE_OUT)
  @test isa(Ωact_out,Gridap.Geometry.BodyFittedTriangulation)
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

  #using Gridap.CellData
  #acell_to_v = get_data(v_in)
  #glue1 = get_glue(Ωact_in,Val(2))
  #cell_to_v = extend(acell_to_v,glue1.mface_to_tface)
  #glue2 = get_glue(Ω_in,Val(2))
  #z0 = cell_to_v.maps.value.values_pos.args[1]
  #print_op_tree(z0)

  #z1 = lazy_map(Reindex(get_array(glue1.mface_to_tface)),glue2.tface_to_mface)
  #z2 = lazy_map(Reindex(z0),z1)

  #print_op_tree(z1)
  #print_op_tree(z2)
  #kkkkk

  #r = lazy_map(Reindex(cell_to_v),glue2.tface_to_mface)
  #print_op_tree(r)
  #kkk

  #rf = change_domain(v_in,Ω_in,ReferenceDomain())
  #print_op_tree(get_data(rf))
  #kkkk

  #r = (∫(v_in)dΩ_in)[Ω_in]
  #print_op_tree(r.args[1])

  #kkk

  # Check divergence theorem
  a = sum( ∫( ∇(v_in)⋅∇(u_in) )*dΩ_in )
  b = sum( ∫( v_in*n_Γ⋅∇(u_in) )*dΓ )
  @test abs(a-b) < 1.0e-9

  scell_val = (∫( ∇(v_in)⋅∇(u_in) )*dΩ_in)[Ω_in]
  @test isa(scell_val,AppendedArray)
  acell_val, Ωa = Gridap.Geometry.move_contributions(scell_val,Ω_in)
  @test isa(acell_val,AppendedArray)
  @test isa(Ωa,AppendedTriangulation)
  sface_val = (∫( v_in*n_Γ⋅∇(u_in) )*dΓ)[Γ]
  aface_val, Γa = Gridap.Geometry.move_contributions(sface_val,Γ)

  # Check divergence theorem (after moving contributions)
  a = sum( acell_val )
  b = sum( aface_val )
  @test abs(a-b) < 1.0e-9

  @test sum(sface_val) ≈ sum(aface_val)
  @test sum(scell_val) ≈ sum(acell_val)

  dv = get_fe_basis(V_in)
  du = get_trial_fe_basis(V_in)
  scell_val = (∫( ∇(v_in)⋅∇(u_in) )*dΩ_in)[Ω_in]
  acell_val, Ωa = Gridap.Geometry.move_contributions(scell_val,Ω_in)
  @test sum(scell_val) ≈ sum(acell_val)

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
