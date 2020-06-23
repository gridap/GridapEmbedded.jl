module BimaterialLinElastCutFEM

using Gridap
using Gridap.Arrays: reindex
using Gridap.Geometry: cell_measure
using GridapEmbedded
using LinearAlgebra: tr

function lame_parameters(E,ν)
  λ = (E*ν)/((1+ν)*(1-2*ν))
  μ = E/(2*(1+ν))
  (λ, μ)
end

function main(;n,outputfile=nothing)

  # Background model
  L = VectorValue(3,2,3)
  partition = (L[1]*n,2*L[2]*n,2*L[3]*n)
  pmin = Point(0.,0.,0.)
  pmax = pmin + L
  bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

  # Identify Dirichlet boundaries
  labeling = get_face_labeling(bgmodel)
  add_tag_from_tags!(labeling,"support0",[1,3,5,7,13,15,17,19,25])
  add_tag_from_tags!(labeling,"support1",[2,4,6,8,14,16,18,20,26])

  # Define geometry
  R = 0.4
  δ = 0.3 + R
  pmid = 0.5*(pmax+pmin)

  geo1 = cylinder(R,x0=Point(0.0,pmid[2],pmin[3]+δ))
  geo2 = cylinder(R,x0=Point(0.0,pmid[2],pmax[3]-δ))

  geo3 = union(geo1,geo2,name="steel")
  geo4 = !(geo3,name="concrete")

  # Material parameters 1
  E1 = 4.0
  ν1 = 0.2
  λ1, μ1 = lame_parameters(E1,ν1)
  @law σ1(ε) = λ1*tr(ε)*one(ε) + 2*μ1*ε

  # Material parameters 2
  E2 = 1.0
  ν2 = 0.2
  λ2, μ2 = lame_parameters(E2,ν2)
  @law σ2(ε) = λ2*tr(ε)*one(ε) + 2*μ2*ε

  # Dirichlet values
  u0 = VectorValue(0,0,0)
  u1 = VectorValue(0.0,0.0,-0.001)

  # Cut the background model
  cutgeo = cut(bgmodel,union(geo3,geo4))

  # Setup models
  model1 = DiscreteModel(cutgeo,"steel")
  model2 = DiscreteModel(cutgeo,"concrete")

  # Setup integration meshes
  trian_Ω1 = Triangulation(cutgeo,"steel")
  trian_Ω2 = Triangulation(cutgeo,"concrete")
  trian_Γ = EmbeddedBoundary(cutgeo,"steel","concrete")

  # Setup normal vectors
  n_Γ = get_normal_vector(trian_Γ)

  # Setup cuadratures
  order = 1
  quad_Ω1 = CellQuadrature(trian_Ω1,2*order)
  quad_Ω2 = CellQuadrature(trian_Ω2,2*order)
  quad_Γ = CellQuadrature(trian_Γ,2*order)

  # Setup stabilization parameters
  meas_K1 = cell_measure(trian_Ω1,num_cells(bgmodel))
  meas_K2 = cell_measure(trian_Ω2,num_cells(bgmodel))
  meas_KΓ = cell_measure(trian_Γ,num_cells(bgmodel))
  
  meas_K1_Γ = reindex(meas_K1,trian_Γ)
  meas_K2_Γ = reindex(meas_K2,trian_Γ)
  meas_KΓ_Γ = reindex(meas_KΓ,trian_Γ)

  γ_hat = 2
  κ1 = (E2*meas_K1_Γ) ./ (E2*meas_K1_Γ .+ E1*meas_K2_Γ)
  κ2 = (E1*meas_K2_Γ) ./ (E2*meas_K1_Γ .+ E1*meas_K2_Γ)
  β = (γ_hat*meas_KΓ_Γ) ./ ( meas_K1_Γ/E1 .+ meas_K2_Γ/E2 )

  # Jump and mean operators for this formulation
  function jump_u(u)
    u1,u2 = u
    u1 - u2
  end
  
  function mean_t(u)
    u1,u2 = u
    κ1*σ1(ε(u1))⋅n_Γ + κ2*σ2(ε(u2))⋅n_Γ
  end
  
  # Setup FESpace
  V1 = TestFESpace(
    model=model1, valuetype=VectorValue{3,Float64},
    reffe=:Lagrangian, order=order, conformity=:H1,
    dirichlet_tags=["support0","support1"])
  
  V2 = TestFESpace(
    model=model2, valuetype=VectorValue{3,Float64},
    reffe=:Lagrangian, order=order, conformity=:H1,
    dirichlet_tags=["support0","support1"])

  U1= TrialFESpace(V1,[u0,u1])
  U2 = TrialFESpace(V2,[u0,u1])

  V = MultiFieldFESpace([V1,V2])
  U = MultiFieldFESpace([U1,U2])

  # Weak form
  function a_Ω1(u,v)
    u1,u2 = u
    v1,v2 = v
    ε(v1) ⊙ σ1(ε(u1))
  end
  
  function a_Ω2(u,v)
    u1,u2 = u
    v1,v2 = v
    ε(v2) ⊙ σ2(ε(u2))
  end

  function a_Γ(u,v)
    β*jump_u(v)⋅jump_u(u) - jump_u(v)⋅mean_t(u) - mean_t(v)⋅jump_u(u)
  end

  # FE problem
  t_Ω1 = LinearFETerm(a_Ω1,trian_Ω1,quad_Ω1)
  t_Ω2 = LinearFETerm(a_Ω2,trian_Ω2,quad_Ω2)
  t_Γ = LinearFETerm(a_Γ,trian_Γ,quad_Γ)
  op = AffineFEOperator(U,V,t_Ω1,t_Ω2,t_Γ)
  uh1, uh2 = solve(op)

  # Postprocess
  uh_Ω1 = restrict(uh1,trian_Ω1)
  uh_Ω2 = restrict(uh2,trian_Ω2)
  if outputfile !== nothing
    writevtk(trian_Ω1,"$(outputfile)_steel",cellfields=["uh"=>uh_Ω1,"sigma"=>σ1(ε(uh_Ω1))])
    writevtk(trian_Ω2,"$(outputfile)_concrete",cellfields=["uh"=>uh_Ω2,"sigma"=>σ2(ε(uh_Ω2))])
  end

end

#main(n=4,outputfile="results")
#main(n=8,outputfile="results")

end # module
