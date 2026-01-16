
using Gridap
using PartitionedArrays
using GridapDistributed

using Gridap.Arrays
using Gridap.Geometry
using Gridap.FESpaces

sol(x) = abs(x[1])
#sol(x) = cos(pi*x[1]/2)

function l2_error(u1,u2)
  Ω = get_triangulation(u1)
  dΩ = Measure(Ω,4)
  eh = u1 - u2
  return sqrt(sum(∫(eh⋅eh)*dΩ))
end

function solve_mass(U,V)
  Ω = get_triangulation(V)
  dΩ = Measure(Ω,4)
  a(u,v) = ∫(u⋅v)*dΩ
  l(v) = ∫(sol*v)*dΩ
  op = AffineFEOperator(a,l,U,V)
  uh = solve(op)
  if isa(uh,GridapDistributed.DistributedCellField)
    consistent!(get_free_dof_values(uh)) |> wait
  end
  eh = uh - sol
  return sqrt(sum(∫(eh⋅eh)*dΩ))
end

function get_solution(U,V)
  Ω = get_triangulation(V)
  dΩ = Measure(Ω,4)
  a(u,v) = ∫(u⋅v)*dΩ
  l(v) = ∫(sol*v)*dΩ
  op = AffineFEOperator(a,l,U,V)
  uh = solve(op)
  if isa(uh,GridapDistributed.DistributedCellField)
    consistent!(get_free_dof_values(uh)) |> wait
  end
  return uh
end

# Serial test

model_serial = CartesianDiscreteModel((-1,1,0,1),(4,1))
V_serial = FESpace(model_serial,ReferenceFE(lagrangian,Float64,1))
Vc_serial = FESpaces.FESpaceWithLinearConstraints(
  Int32[5,10],Table([[1],[6]]),Table([[1.0],[1.0]]),V_serial
)
solve_mass(Vc_serial,Vc_serial)
solve_mass(V_serial,V_serial)

# Distributed test

np = (2,1)
ranks = LinearIndices((prod(np),))

model = CartesianDiscreteModel(ranks,np,(-1,1,0,1),(4,1))

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(model, reffe)
solve_mass(V,V)

############################################################################################

dof_lids = get_cell_dof_ids(V)
dof_gids = partition(get_free_dof_ids(V))

spaces, ids = map(ranks,local_views(V)) do r, V
  if r == 1
    n_fmdofs = 8
    mDOF_to_dof = Int32[2,3,4,6,7,8,0,0]
    sDOF_to_dof = Int32[1,5]
    sDOF_to_mdofs = Table([[7],[8]])
    sDOF_to_coeffs = Table([[1.0],[1.0]])

    ids = LocalIndices(
      8, r, Int[1,3,4,2,6,7,5,8], Int32[1,2,2,1,2,2,2,2]
    )
  else
    n_fmdofs = num_free_dofs(V)
    mDOF_to_dof = collect(Int32,1:n_fmdofs)
    sDOF_to_dof = Int32[]
    sDOF_to_mdofs = empty_table(Int32,Int,0)
    sDOF_to_coeffs = empty_table(Float64,Int,0)
  
    ids = LocalIndices(
      8, r, Int[1,3,4,5,2,6,7,8], Int32[1,2,2,2,1,2,2,2]
    )
  end
  Vc = FESpaces.FESpaceWithLinearConstraints(
    V, mDOF_to_dof, sDOF_to_dof, sDOF_to_mdofs, sDOF_to_coeffs, n_fmdofs
  )
  return Vc, ids
end |> tuple_of_arrays

gids = PRange(ids)
trian = GridapDistributed.DistributedTriangulation(map(get_triangulation,spaces),model)
vector_type = GridapDistributed._find_vector_type(spaces,gids)
Vc = GridapDistributed.DistributedSingleFieldFESpace(spaces,gids,trian,vector_type)

solve_mass(Vc,Vc)


u1 = interpolate(sol,V)
u2 = interpolate(sol,Vc)
x1 = get_free_dof_values(u1)
x2 = get_free_dof_values(u2)
consistent!(x2) |> wait

l2_error(u1,u2)

uhc = get_solution(Vc,Vc)
writevtk(Triangulation(model), "data/dist", cellfields=["u"=>sol,"uh"=>uhc, "eh"=>(uhc - sol)])
