
using Gridap
using GridapEmbedded
using Gridap.Geometry, Gridap.Arrays
using Gridap.FESpaces, Gridap.CellData

function get_aggregates(cutgeo)
  bgmodel = get_background_model(cutgeo)
  geo = get_geometry(cutgeo)
  strategy = AggregateAllCutCells()
  aggregates = aggregate(strategy,cutgeo,geo)
  colors = color_aggregates(aggregates,bgmodel)
  return inverse_table(colors)
end

function projection_operator(ptopo, V, W, Ω, dΩ)
  Π(u,Ω) = change_domain(u,Ω,DomainStyle(u))
  mass(u,v) = ∫(u⋅Π(v,Ω))dΩ
  P = LocalOperator(
    LocalSolveMap(), ptopo, W, V, mass, mass
  )
  return P
end

R = 1.4
geo = disk(R;x0=Point(1.0,1.7))

n = 8
order = 1
qdegree = 2*order

bgmodel = CartesianDiscreteModel((0,1,0,1),(n,n))
cutgeo = cut(bgmodel,geo)

Ω_ac = Triangulation(cutgeo,ACTIVE)
V = FESpace(Ω_ac,ReferenceFE(raviart_thomas,Float64,order))
Q = FESpace(Ω_ac,ReferenceFE(lagrangian,Float64,order); conformity=:L2)
X = MultiFieldFESpace([V,Q])

ptopo = PatchTopology(get_grid_topology(bgmodel), get_aggregates(cutgeo))
Ωp = PatchTriangulation(bgmodel, ptopo)
dΩp = Measure(Ωp, qdegree)

Wd = FESpaces.PatchFESpace(bgmodel,Ωp,VectorValue{2,Float64},order)
Πd = projection_operator(ptopo, V, W, Ωp, dΩp)
W0 = FESpaces.PatchFESpace(bgmodel,Ωp,Float64,order)
Π0 = projection_operator(ptopo, Q, W0, Ωp, dΩp)

Ω = Triangulation(cutgeo,PHYSICAL)
Ω_cut = Triangulation(cutgeo,CUT_IN)
Γ = EmbeddedBoundary(cutgeo)
dΩ = Measure(Ω, qdegree)
dΩ_cut = Measure(Ω_cut, qdegree)
dΓ = Measure(Γ, 2*qdegree)

η = 1.0
u_exact(x) = VectorValue(x[1], -x[2])
p_exact(x) = x[1]*x[2]
f(x) = η*u_exact(x) + ∇(p_exact)(x)
g(x) = -(∇⋅u_exact)(x)

γ = 10.0
hΓ = CellField(γ ./ get_array(∫(1.0)dΓ),Γ)
n_Γ = get_normal_vector(Γ)
a(u,v) = ∫(u⋅v)dΩ + ∫(hΓ*(u⋅n_Γ)*(v⋅n_Γ))dΓ
b(u,q) = ∫((∇⋅u)*q)dΩ
btilde(u,q) = b(u,q) + ∫((u⋅n_Γ)*q)dΓ
l((v,q)) = ∫(v⋅f + g⋅q)dΩ - ∫((v⋅n_Γ)*p_exact)dΓ + ∫(hΓ*(v⋅n_Γ)⋅(u_exact⋅n_Γ))dΓ

τ = 1.0
sd(u,v) = ∫(τ*(u⋅v))dΩ_cut
s0(p,q) = ∫(τ*p*q)dΩ_cut

function weakform(x,y)
  u, p = x
  v, q = y
  Πu, Πv = Π(u), Π(v)
  Xp = FESpaces.PatchFESpace(X,ptopo)
  data = FESpaces.collect_and_merge_cell_matrix_and_vector(
    (X, X  , a(u,v) + btilde(v,p) + b(u,q) + s(u,v), l(y), zero(X)),
    (X, Xp , -1*s(u,Πv), DomainContribution(), zero(X)),
    (Xp, X , -1*s(Πu,v), DomainContribution(), zero(Xp)),
    (Xp, Xp, s(Πu,Πv), DomainContribution(), zero(Xp)),
  )
  assem = SparseMatrixAssembler(X,X)
  A, B = assemble_matrix_and_vector(assem,data)
  return AffineFEOperator(X,X,A,B)
end

op = weakform(get_trial_fe_basis(X),get_fe_basis(X))
uh, ph = solve(op)

l2_err_u = sqrt(sum(∫( (uh-u_exact)⋅(uh-u_exact) )dΩ))
l2_err_p = sqrt(sum(∫( (ph-p_exact)*(ph-p_exact) )dΩ))
