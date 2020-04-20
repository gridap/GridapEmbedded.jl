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
  trian_Ω = Triangulation(cutgeo,"csg")
  trian_Γd = EmbeddedBoundary(cutgeo,"csg","source")
  trian_Γg = GhostSkeleton(cutgeo,"csg")

  # Setup normal vectors
  n_Γd = get_normal_vector(trian_Γd)
  n_Γg = get_normal_vector(trian_Γg)

  #writevtk(trian_Ω,"trian_O")
  #writevtk(trian_Γd,"trian_Gd",cellfields=["normal"=>n_Γd])
  #writevtk(trian_Γg,"trian_Gg",cellfields=["normal"=>n_Γg])
  #writevtk(Triangulation(bgmodel),"bgtrian")

  # Setup cuadratures
  order = 1
  quad_Ω = CellQuadrature(trian_Ω,2*order)
  quad_Γd = CellQuadrature(trian_Γd,2*order)
  quad_Γg = CellQuadrature(trian_Γg,2*order)

  # Setup FESpace
  model = DiscreteModel(cutgeo,"csg")
  V = TestFESpace(model=model,valuetype=Float64,reffe=:Lagrangian,order=order,conformity=:H1)
  U = TrialFESpace(V)

  # Weak form
  γd = 10.0
  γg = 0.1
  h = (pmax - pmin)[1] / partition[1]
  a_Ω(u,v) = ∇(v)*∇(u)
  l_Ω(v) = v*f
  a_Γd(u,v) = (γd/h)*v*u  - v*(n_Γd*∇(u)) - (n_Γd*∇(v))*u
  l_Γd(v) = (γd/h)*v*ud - (n_Γd*∇(v))*ud
  a_Γg(v,u) = (γg*h)*jump(n_Γg*∇(v))*jump(n_Γg*∇(u))

  # FE problem
  t_Ω = AffineFETerm(a_Ω,l_Ω,trian_Ω,quad_Ω)
  t_Γd = AffineFETerm(a_Γd,l_Γd,trian_Γd,quad_Γd)
  t_Γg = LinearFETerm(a_Γg,trian_Γg,quad_Γg)
  op = AffineFEOperator(U,V,t_Ω,t_Γd,t_Γg)
  uh = solve(op)

  # Postprocess
  uh_Ω = restrict(uh,trian_Ω)
  if outputfile !== nothing
    writevtk(trian_Ω,outputfile,cellfields=["uh"=>uh_Ω])
  end

end

main(n=40,outputfile="results")

end # module
