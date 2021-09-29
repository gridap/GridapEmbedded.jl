module StokesTubeWithObstacleCutFEM

using Gridap
import Gridap: ∇
using GridapEmbedded
using Test
using LinearAlgebra: tr

function main(;n,outputfile=nothing)

  # Formulation taken from
  # André Massing · Mats G. Larson · Anders Logg · Marie E. Rognes,
  # A Stabilized Nitsche Fictitious Domain Method for the Stokes Problem
  # J Sci Comput (2014) 61:604–628 DOI 10.1007/s10915-014-9838-9

  # Select geometry
  R = 0.3
  L = 2.0
  x0 = Point(0.0,0.0,0.0)
  geo2 = tube(R,L,x0=x0)
  geo3 = sphere(0.3*R,x0=x0+Point(L/3,0.0,0.0))
  geo4 = ! get_geometry(geo2,"walls")
  geo5 = union(geo3,geo4,name="solid")
  geo1 = setdiff(geo2,geo3,name="fluid")

   m = 0.05
  pmin = x0-Point(m,R,R)
  pmax = x0+Point(m+L,R,R)

  # Forcing data
  function uin(x,x0,R,vmax)
    dx = x-x0
    r2 = dx⋅dx
    if r2 < R^2
      v =  (1-r2/R^2)*vmax
    else
      v = zero(vmax)
    end
    v
  end

  vmax = VectorValue(1.0,0.0,0.0)
  uin(x) = uin(x,x0,R,vmax)

  # Background model
  partition = (4*n,n,n)
  D = length(partition)
  bgmodel = simplexify(CartesianDiscreteModel(pmin,pmax,partition))

  # Cut the background model
  cutgeo = cut(bgmodel,union(geo1,geo5))

  # Generate the "active" mesh
  Ω_act = Triangulation(cutgeo,ACTIVE,"fluid")

  # Setup integration meshes
  Ω = Triangulation(cutgeo,PHYSICAL,"fluid")
  Γi = EmbeddedBoundary(cutgeo,"fluid","inlet")
  Γw = EmbeddedBoundary(cutgeo,"fluid","solid")
  Γg = GhostSkeleton(cutgeo,"fluid")

  # Setup normal vectors
  n_Γi = get_normal_vector(Γi)
  n_Γw = get_normal_vector(Γw)
  n_Γg = get_normal_vector(Γg)

  #writevtk(Ω,"trian_O")
  #writevtk(Γi,"trian_Gi",cellfields=["uin"=>uin,"normal"=>n_Γi])
  #writevtk(Γw,"trian_Gw",cellfields=["normal"=>n_Γw])
  #writevtk(Triangulation(bgmodel),"bgtrian")

  # Setup Lebesgue measures
  order = 1
  degree = 2*order
  dΩ = Measure(Ω,degree)
  dΓi = Measure(Γi,degree)
  dΓw = Measure(Γw,degree)
  dΓg = Measure(Γg,degree)

  # Setup FESpace

  reffe_u = ReferenceFE(lagrangian,VectorValue{D,Float64},order,space=:P)
  reffe_p = ReferenceFE(lagrangian,Float64,order,space=:P)

  V = TestFESpace(Ω_act,reffe_u,conformity=:H1)
  Q = TestFESpace(Ω_act,reffe_p,conformity=:H1)

  U = TrialFESpace(V)
  P = TrialFESpace(Q)

  X = MultiFieldFESpace([U,P])
  Y = MultiFieldFESpace([V,Q])

  # Stabilization parameters
  β0 = 0.25
  β1 = 0.2
  β2 = 0.1
  β3 = 0.05
  γ = 10.0
  h = (pmax-pmin)[1]/partition[1]

  # Weak form
  a_Ω(u,v) = ∇(u)⊙∇(v)
  b_Ω(v,p) = - (∇⋅v)*p
  c_Ω(p,q) = (β1*h^2)*∇(p)⋅∇(q)
  a_Γ(u,v,n_Γ) = - (n_Γ⋅∇(u))⋅v - u⋅(n_Γ⋅∇(v)) + (γ/h)*u⋅v
  b_Γ(v,p,n_Γ) = (n_Γ⋅v)*p
  i_Γg(u,v) = (β2*h)*jump(n_Γg⋅∇(u))⋅jump(n_Γg⋅∇(v))
  j_Γg(p,q) = (β3*h^3)*jump(n_Γg⋅∇(p))*jump(n_Γg⋅∇(q))

  a((u,p),(v,q)) =
    ∫( a_Ω(u,v)+b_Ω(u,q)+b_Ω(v,p)-c_Ω(p,q) ) * dΩ +
    ∫( a_Γ(u,v,n_Γi)+b_Γ(u,q,n_Γi)+b_Γ(v,p,n_Γi) ) * dΓi +
    ∫( a_Γ(u,v,n_Γw)+b_Γ(u,q,n_Γw)+b_Γ(v,p,n_Γw) ) * dΓw +
    ∫( i_Γg(u,v) - j_Γg(p,q) ) * dΓg

  l((v,q)) = ∫( uin⊙( (γ/h)*v - n_Γi⋅∇(v) + q*n_Γi ) ) * dΓi

  # FE problem
  op = AffineFEOperator(a,l,X,Y)

  uh, ph = solve(op)

  # Postprocess
  if outputfile !== nothing
    writevtk(Ω,outputfile, cellfields=["uh"=>uh,"ph"=>ph])
  end

end

end # module
