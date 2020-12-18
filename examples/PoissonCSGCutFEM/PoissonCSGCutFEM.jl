module PoissonCSGCutFEM

using Gridap
using GridapEmbedded

function main(;n,outputfile=nothing)

  # Background model
  partition = (n,n,n)
  pmin = 0.8*Point(-1,-1,-1)
  pmax = 0.8*Point(1,1,1)
  bgmodel = CartesianDiscreteModel(pmin,pmax,partition)

  # Select geometry
  R = 0.5
  geo1 = cylinder(R,v=VectorValue(1,0,0))
  geo2 = cylinder(R,v=VectorValue(0,1,0))
  geo3 = cylinder(R,v=VectorValue(0,0,1))
  geo4 = union(union(geo1,geo2),geo3,name="source")
  geo5 = sphere(1)
  geo6 = cube(L=1.5)
  geo7 = intersect(geo6,geo5)
  geo8 = setdiff(geo7,geo4,name="csg")

  # Forcing data
  ud = 1
  f = 10

  # Cut the background model
  cutgeo = cut(bgmodel,geo8)

  # Setup integration meshes
  Ω = Triangulation(cutgeo,"csg")
  Γd = EmbeddedBoundary(cutgeo,"csg","source")
  Γg = GhostSkeleton(cutgeo,"csg")

  # Setup normal vectors
  n_Γd = get_normal_vector(Γd)
  n_Γg = get_normal_vector(Γg)

  #writevtk(Ω,"trian_O")
  #writevtk(Γd,"trian_Gd",cellfields=["normal"=>n_Γd])
  #writevtk(Γg,"trian_Gg",cellfields=["normal"=>n_Γg])
  #writevtk(Triangulation(bgmodel),"bgtrian")

  # Setup Lebesgue measures
  order = 1
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓd = Measure(Γd,degree)
  dΓg = Measure(Γg,degree)

  # Setup FESpace
  model = DiscreteModel(cutgeo,"csg")
  V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order),conformity=:H1)
  U = TrialFESpace(V)

  # Weak form
  γd = 10.0
  γg = 0.1
  h = (pmax - pmin)[1] / partition[1]

  a(u,v) =
    ∫( ∇(v)⋅∇(u) ) * dΩ +
    ∫( (γd/h)*v*u  - v*(n_Γd⋅∇(u)) - (n_Γd⋅∇(v))*u ) * dΓd +
    ∫( (γg*h)*jump(n_Γg⋅∇(v))*jump(n_Γg⋅∇(u)) ) * dΓg

  l(v) =
    ∫( v*f ) * dΩ +
    ∫( (γd/h)*v*ud - (n_Γd⋅∇(v))*ud ) * dΓd

  # FE problem
  op = AffineFEOperator(a,l,U,V)
  uh = solve(op)

  # Postprocess
  if outputfile !== nothing
    writevtk(Ω,outputfile,cellfields=["uh"=>uh])
  end

end

end # module
