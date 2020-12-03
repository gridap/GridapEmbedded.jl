module BimaterialLinElastCutFEM

using Gridap
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
  σ1(ε) = λ1*tr(ε)*one(ε) + 2*μ1*ε

  # Material parameters 2
  E2 = 1.0
  ν2 = 0.2
  λ2, μ2 = lame_parameters(E2,ν2)
  σ2(ε) = λ2*tr(ε)*one(ε) + 2*μ2*ε

  # Dirichlet values
  u0 = VectorValue(0,0,0)
  u1 = VectorValue(0.0,0.0,-0.001)

  # Cut the background model
  cutgeo = cut(bgmodel,union(geo3,geo4))

  # Setup models
  model1 = DiscreteModel(cutgeo,"steel")
  model2 = DiscreteModel(cutgeo,"concrete")

  # Setup integration meshes
  Ω1 = Triangulation(cutgeo,"steel")
  Ω2 = Triangulation(cutgeo,"concrete")
  Γ = EmbeddedBoundary(cutgeo,"steel","concrete")

  # Setup normal vectors
  n_Γ = get_normal_vector(Γ)

  # Setup Lebesgue measures
  order = 1
  degree = 2*order
  dΩ1 = LebesgueMeasure(Ω1,degree)
  dΩ2 = LebesgueMeasure(Ω2,degree)
  dΓ = LebesgueMeasure(Γ,degree)

  # Setup FESpace

  V1 = TestFESpace(model1,
                   ReferenceFE(:Lagrangian,VectorValue{3,Float64},order),
                   conformity=:H1,
                   dirichlet_tags=["support0","support1"])
  V2 = TestFESpace(model2,
                   ReferenceFE(:Lagrangian,VectorValue{3,Float64},order),
                   conformity=:H1,
                   dirichlet_tags=["support0","support1"])

  U1 = TrialFESpace(V1,[u0,u1])
  U2 = TrialFESpace(V2,[u0,u1])

  V = MultiFieldFESpace([V1,V2])
  U = MultiFieldFESpace([U1,U2])

  # Setup stabilization parameters

  meas_K1 = get_cell_measure(Ω1)
  meas_K2 = get_cell_measure(Ω2)
  meas_KΓ = get_cell_measure(Γ)

  γ_hat = 2
  κ1 = (E2*meas_K1) ./ (E2*meas_K1 .+ E1*meas_K2)
  κ2 = (E1*meas_K2) ./ (E2*meas_K1 .+ E1*meas_K2)
  β = (γ_hat*meas_KΓ) ./ ( meas_K1/E1 .+ meas_K2/E2 )

  # Jump and mean operators for this formulation

  jump_u(u1,u2) = u1 - u2
  mean_t(u1,u2) = κ1*(σ1∘ε(u1)) + κ2*(σ2∘ε(u2))

  # Weak form

  a((u1,u2),(v1,v2)) =
    ∫( ε(v1) ⊙ (σ1∘ε(u1)) ) * dΩ1 + ∫( ε(v2) ⊙ (σ2∘ε(u2)) ) * dΩ2 +
    ∫( β*jump_u(v1,v2)⋅jump_u(u1,u2)
       - n_Γ⋅mean_t(u1,u2)⋅jump_u(v1,v2)
       - n_Γ⋅mean_t(v1,v2)⋅jump_u(u1,u2) ) * dΓ

  l((v1,v2)) = 0

  # FE problem
  op = AffineFEOperator(a,l,U,V)
  uh1, uh2 = solve(op)
  uh = (uh1,uh2)

  # Postprocess
  if outputfile !== nothing
    writevtk(Ω1,"$(outputfile)_steel",cellfields=["uh"=>uh1,"sigma"=>σ1∘ε(uh1)])
    writevtk(Ω2,"$(outputfile)_concrete",cellfields=["uh"=>uh2,"sigma"=>σ2∘ε(uh2)])
  end

end

#main(n=4,outputfile="results")
#main(n=8,outputfile="results")

end # module
