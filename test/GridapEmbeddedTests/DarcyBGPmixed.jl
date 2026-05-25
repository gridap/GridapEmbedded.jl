module DarcyBGPmixed
# Darcy problem: BGP stabilisation with mixed boundary conditions (flux + pressure traces
# imposed on different parts of the unfitted boundary).
#
# Ref: Badia, Boschman, Martín, Nilsson, Ruiz-Baier, Zahedi (2026)
#      "Divergence-free unfitted finite element discretisations for the Darcy problem"
#      arXiv:2603.26212

using Test
using Gridap
using GridapEmbedded
using Gridap.Geometry, Gridap.Arrays
using Gridap.FESpaces, Gridap.CellData

function get_aggregates(cutgeo)
  geo = get_geometry(cutgeo)
  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo,geo)
  return Table(filter(agg -> length(agg) > 1, inverse_table(aggregates)))
end

function projection_operator(ptopo, V, W, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  P = LocalOperator(
    LocalSolveMap(), ptopo, W, V, mass, mass
  )
  return P
end

function div_projection_operator(ptopo, V, W, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  lhs(u,v) = ∫(u⋅Π(v,Ω))dΩ
  rhs(u,v) = ∫((∇⋅u)⋅Π(v,Ω))dΩ
  P = LocalOperator(
    LocalSolveMap(), ptopo, W, V, lhs, rhs
  )
  return P
end

function generate_mesh(n,R=0.4,x0=0.835)
  geo1 = disk(R;x0=Point(0.5,0.5),name="disk")
  geo2 = plane(x0=Point(x0,0.5),v=Point(1.0,0.0),name="plane")
  geo = intersect(geo1,geo2,name="domain")

  bgmodel = CartesianDiscreteModel((0,1,0,1),(n,n))
  cutgeo = cut(bgmodel,geo)
  return cutgeo
end

function driver(n,order,u_exact,p_exact,η=1.0,γ=10.0,τ=10.0)
  cutgeo = generate_mesh(n)
  bgmodel = get_background_model(cutgeo)
  qdegree = 2*(order+1)

  Ω_ac = Triangulation(cutgeo,ACTIVE)
  V = FESpace(Ω_ac,ReferenceFE(raviart_thomas,Float64,order); conformity=:Hdiv)
  Q = FESpace(Ω_ac,ReferenceFE(lagrangian,Float64,order); conformity=:L2)
  X = MultiFieldFESpace([V,Q])

  ptopo = PatchTopology(get_grid_topology(bgmodel), get_aggregates(cutgeo))
  Ωp = PatchTriangulation(bgmodel, ptopo)
  
  Ω = Triangulation(cutgeo,PHYSICAL)
  Γu = EmbeddedBoundary(cutgeo, "disk","domain")
  Γp = EmbeddedBoundary(cutgeo, "plane","domain")

  dΩ = Measure(Ω, qdegree)
  dΩp = Measure(Ωp, qdegree)
  dΓu = Measure(Γu, 2*qdegree)
  dΓp = Measure(Γp, 2*qdegree)

  Wd = FESpaces.PatchFESpace(bgmodel,Ωp,Float64,order;space=:RT)
  Πd = projection_operator(ptopo, V, Wd, Ωp, dΩp)
  W0 = FESpaces.PatchFESpace(bgmodel,Ωp,Float64,order;space=:Q)
  Π0 = projection_operator(ptopo, Q, W0, Ωp, dΩp)
  Π0div = div_projection_operator(ptopo, V, W0, Ωp, dΩp)

  f(x) = η*u_exact(x) + ∇(p_exact)(x)
  g(x) = -(∇⋅u_exact)(x)

  h = γ/n
  n_Γu = get_normal_vector(Γu)
  n_Γp = get_normal_vector(Γp)
  a(u,v) = ∫(η*(u⋅v))dΩ + ∫(h*(u⋅n_Γu)*(v⋅n_Γu))dΓu
  b(u,q) = ∫(-(∇⋅u)*q)dΩ
  btilde(u,q) = b(u,q) + ∫((u⋅n_Γu)*q)dΓu
  l((v,q)) = ∫(v⋅f + g⋅q)dΩ  + ∫(h*(v⋅n_Γu)⋅(u_exact⋅n_Γu))dΓu - ∫((v⋅n_Γp)*p_exact)dΓp

  sd(u,v) = ∫(τ*(u⋅v))dΩp
  s0(p,q) = ∫(τ*p*q)dΩp

  function weakform(x,y)
    u, p = x
    v, q = y
    Πu, Πv = Πd(u), Πd(v)
    Πp, Πq = Π0(p), Π0(q)
    divu, divv = ∇⋅u, ∇⋅v
    Πdivu, Πdivv = Π0div(u), Π0div(v)
    Xp = FESpaces.PatchFESpace(X,ptopo)
    data = FESpaces.collect_and_merge_cell_matrix_and_vector(
      (X, X  , a(u,v) + btilde(v,p) + b(u,q) + sd(u,v) - s0(divu,q) - s0(divv,p), l(y)),
      (X, Xp , s0(divu,Πq) + s0(Πdivv,p) - sd(u,Πv), DomainContribution()),
      (Xp, X , s0(Πdivu,q) + s0(divv,Πp) - sd(Πu,v), DomainContribution()),
      (Xp, Xp, sd(Πu,Πv) - s0(Πdivu,Πq) - s0(Πdivv,Πp), DomainContribution()),
    )
    assem = SparseMatrixAssembler(X,X)
    A, B = assemble_matrix_and_vector(assem,data)
    return AffineFEOperator(X,X,A,B)
  end

  op = weakform(get_trial_fe_basis(X),get_fe_basis(X))
  uh, ph = solve(op)

  eu = uh-u_exact
  ep = ph-p_exact
  l2_err_u = sqrt(sum(∫( eu ⋅ eu )dΩ))
  l2_err_p = sqrt(sum(∫( ep ⋅ ep )dΩ))

  return cutgeo, uh, ph, l2_err_u, l2_err_p
end

function convergence(ns,order,u_exact,p_exact)
  l2_u = zeros(Float64,length(ns))
  l2_p = zeros(Float64,length(ns))
  for (i,n) in enumerate(ns)
    _, _, _, l2_u[i], l2_p[i] = driver(n,order,u_exact,p_exact)
  end
  slope_u = log.(l2_u[1:end-1]./l2_u[2:end]) ./ log.(ns[2:end]./ns[1:end-1])
  slope_p = log.(l2_p[1:end-1]./l2_p[2:end]) ./ log.(ns[2:end]./ns[1:end-1])
  return l2_u, l2_p, slope_u, slope_p
end

# Manufactured solution
u_exact(x) = VectorValue(x[1], -x[2])
p_exact(x) = x[1] + x[2]
cutgeo, uh, ph, l2_err_u, l2_err_p = driver(8,1,u_exact,p_exact)
@test l2_err_u < 1.e-10
@test l2_err_p < 1.e-10

u_conv(x) = VectorValue(x[1] + sin(π*x[2]), -x[2] + sin(π*x[1]))
p_conv(x) = sin(π*x[1]) - sin(π*x[2])
# l2_u, l2_p, su, sp = convergence([8,16,32,64],1,u_conv,p_conv)

writevtk(
  Ω, "darcy_BGP_mixed"; append=false, 
  cellfields=[
    "uh"=>uh,"ph"=>ph,"u_exact"=>u_exact,"p_exact"=>p_exact,
    "eu"=>uh-u_exact,"ep"=>ph-p_exact
  ],
)

cell_to_agg = flatten_partition(get_aggregates(cutgeo),num_cells(bgmodel))
writevtk(
  Triangulation(bgmodel), "aggregates"; append=false,
  celldata = ["agg"=>cell_to_agg]
)

function normal_at_centroid(Γ)
  n = get_normal_vector(Γ)
  x = CellPoint(map(mean,get_cell_coordinates(Γ)),Γ,PhysicalDomain())
  return n(x)
end

Γu = EmbeddedBoundary(cutgeo, "disk","domain")
Γp = EmbeddedBoundary(cutgeo, "plane","domain")
writevtk(
  Γu, "boundary_u"; append=false, celldata=["n"=>normal_at_centroid(Γu)]
)
writevtk(
  Γp, "boundary_p"; append=false, celldata=["n"=>normal_at_centroid(Γp)]
)

end # module