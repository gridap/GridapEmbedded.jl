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
    r2 = dx*dx
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
  n = 10
  partition = (4*n,n,n)
  D = length(partition)
  bgmodel = simplexify(CartesianDiscreteModel(pmin,pmax,partition))
  
  # Cut the background model
  cutgeo = cut(bgmodel,union(geo1,geo5))
  
  # Generate the "active" model
  model = DiscreteModel(cutgeo,"fluid")
  
  # Setup integration meshes
  trian_Ω = Triangulation(cutgeo,"fluid")
  trian_Γi = EmbeddedBoundary(cutgeo,"fluid","inlet")
  trian_Γw = EmbeddedBoundary(cutgeo,"fluid","solid")
  trian_Γg = GhostSkeleton(cutgeo,"fluid")
  
  # Setup normal vectors
  n_Γi = get_normal_vector(trian_Γi)
  n_Γw = get_normal_vector(trian_Γw)
  n_Γg = get_normal_vector(trian_Γg)
  
  #writevtk(trian_Ω,"trian_O")
  #writevtk(trian_Γi,"trian_Gi",cellfields=["uin"=>uin,"normal"=>n_Γi])
  #writevtk(trian_Γw,"trian_Gw",cellfields=["normal"=>n_Γw])
  #writevtk(Triangulation(bgmodel),"bgtrian")
  
  # Setup cuadratures
  order = 1
  quad_Ω = CellQuadrature(trian_Ω,2*order)
  quad_Γi = CellQuadrature(trian_Γi,2*order)
  quad_Γw = CellQuadrature(trian_Γw,2*order)
  quad_Γg = CellQuadrature(trian_Γg,2*order)
  
  # Setup FESpace
  
  V = TestFESpace(
    model=model,valuetype=VectorValue{D,Float64},reffe=:PLagrangian,
    order=order,conformity=:H1)
  
  Q = TestFESpace(
    model=model,valuetype=Float64,reffe=:PLagrangian,
    order=order,conformity=:H1)
  
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
  a_Ω(u,v) = inner(∇(u),∇(v))
  b_Ω(v,p) = - (∇*v)*p
  c_Ω(p,q) = (β1*h^2)*∇(p)*∇(q)
  a_Γ(u,v,n_Γ) = - (n_Γ*∇(u))*v - u*(n_Γ*∇(v)) + (γ/h)*u*v
  b_Γ(v,p,n_Γ) = (n_Γ*v)*p
  i_Γg(u,v) = (β2*h)*jump(n_Γg*∇(u))*jump(n_Γg*∇(v))
  j_Γg(p,q) = (β3*h^3)*jump(n_Γg*∇(p))*jump(n_Γg*∇(q))
  
  function A_Ω(X,Y)
    u,p = X
    v,q = Y
    a_Ω(u,v)+b_Ω(u,q)+b_Ω(v,p)-c_Ω(p,q)
  end
  
  function A_Γ(X,Y,n_Γ)
    u,p = X
    v,q = Y
    a_Γ(u,v,n_Γ)+b_Γ(u,q,n_Γ)+b_Γ(v,p,n_Γ)
  end
  
  function J_Γg(X,Y)
    u,p = X
    v,q = Y
    i_Γg(u,v) - j_Γg(p,q) 
  end
  
  function L_Γi(Y)
    v,q = Y
    uin*( (γ/h)*v - n_Γi*∇(v) + q*n_Γi )
  end
  
  # FE problem
  t_Ω = LinearFETerm(A_Ω,trian_Ω,quad_Ω)
  t_Γi = AffineFETerm((X,Y)->A_Γ(X,Y,n_Γi),L_Γi,trian_Γi,quad_Γi)
  t_Γw = LinearFETerm((X,Y)->A_Γ(X,Y,n_Γw),trian_Γw,quad_Γw)
  t_Γg = LinearFETerm(J_Γg,trian_Γg,quad_Γg)
  op = AffineFEOperator(X,Y,t_Ω,t_Γi,t_Γw,t_Γg)
  uh, ph = solve(op)
  
  # Postprocess
  uh_Ω = restrict(uh,trian_Ω)
  ph_Ω = restrict(ph,trian_Ω)
  if outputfile !== nothing
    writevtk(trian_Ω,outputfile, cellfields=["uh"=>uh_Ω,"ph"=>ph_Ω])
  end
  
end

end # module

