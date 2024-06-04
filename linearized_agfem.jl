using Gridap
using GridapEmbedded

function u_ex(x)
   x[1]+x[2]
end

f(x)=-Δ(u_ex)(x)

function Gridap.Geometry.get_background_model(
     a::GridapEmbedded.Interfaces.SubCellTriangulation{Dc,
                                                       Dp,
                                                       T,
                                                       <:Gridap.Adaptivity.AdaptedDiscreteModel}) where {Dc,Dp,T}
  Gridap.Adaptivity.get_model(a.bgmodel)
end

const SFC_ADM{Dc,Dp,T} =  
     GridapEmbedded.Interfaces.SubFacetTriangulation{Dc,Dp,T,<:Gridap.Adaptivity.AdaptedDiscreteModel}

function Gridap.Geometry.get_background_model(
  a::SFC_ADM{Dc,Dp,T,<:Gridap.Adaptivity.AdaptedDiscreteModel}) where {Dc,Dp,T}
  Gridap.Adaptivity.get_model(a.bgmodel)
end

const APP_ST_AT{Dc,Dp} =  Gridap.Geometry.AppendedTriangulation{Dc,Dp,
                                                               <:GridapEmbedded.Interfaces.SubCellTriangulation,
                                                               <:Gridap.Adaptivity.AdaptedTriangulation}

const APP_AT_AT{Dc,Dp} =  Gridap.Geometry.AppendedTriangulation{Dc,Dp,
                                                               <:Gridap.Adaptivity.AdaptedTriangulation,
                                                               <:Gridap.Adaptivity.AdaptedTriangulation}

function Gridap.FESpaces.get_cell_fe_data(fun,f,ttrian::Union{APP_ST_AT{Dc,Dp},APP_AT_AT{Dc,Dp}}) where {Dc,Dp}
  sface_to_data = fun(f)
  strian = get_triangulation(f)
  Gridap._get_cell_fe_data(fun,sface_to_data, strian, ttrian)
end

function Gridap.FESpaces.get_cell_fe_data(fun,f,ttrian::SFC_ADM{Dc,Dp,T}) where {Dc,Dp,T}
  sface_to_data = fun(f)
  strian = get_triangulation(f)
  Gridap._get_cell_fe_data(fun,sface_to_data, strian, ttrian)
end

function  Gridap._get_cell_fe_data(fun, 
  sface_to_data, 
  strian::Gridap.Geometry.BodyFittedTriangulation{Dc,Dp}, 
  ttrian::Union{APP_ST_AT{Dc,Dp},APP_AT_AT{Dc,Dp}}) where {Dc,Dp}
  
  bt_adapted_model = Gridap.Adaptivity.get_adapted_model(ttrian.b)
  if (get_background_model(strian) === Gridap.Adaptivity.get_parent(bt_adapted_model))
  #  # Fine-grid background model cells that are in+cut cells 
  #  tglue=get_glue(ttrian,Val(Dc))
  #  m_fine_cells_in_cut = tglue.tface_to_mface
  #  # Coarse-grid background models cells associated to fine-grid cells that are in+cut cells 
  #  m_coarse_cells_in_cut = bt_adapted_model.glue.n2o_faces_map[end][m_fine_cells_in_cut]
  #   # Coarse-grid triangulation cells that are in+cut cells
  #  sglue=get_glue(strian,Val(Dc))
  #  t_coarse_cells_in_cut = sglue.mface_to_tface[m_coarse_cells_in_cut]
  #  lazy_map(Reindex(sface_to_data),t_coarse_cells_in_cut)
    tglue = get_glue(ttrian,Val(Dc))
    ttrian_tface_to_mface=tglue.tface_to_mface
    strian_ttrian_tface_to_mface=bt_adapted_model.glue.n2o_faces_map[end][ttrian_tface_to_mface]
    sglue = get_glue(strian,Val(Dc))
    strian_mface_to_tface=sglue.mface_to_tface
    lazy_map(Reindex(sface_to_data),strian_mface_to_tface[strian_ttrian_tface_to_mface])
  elseif (get_background_model(strian) === get_background_model(bt))
    Gridap.Helpers.@notimplemented
  else
    Gridap.Helpers.@notimplemented
  end
end

function Gridap._get_cell_fe_data(fun, 
  sface_to_data, 
  strian::Gridap.Geometry.BodyFittedTriangulation{Dc,Dp}, 
  ttrian::SFC_ADM{Df,Dp}) where {Dc,Df,Dp}
  
  if (get_background_model(strian) === Gridap.Adaptivity.get_parent(ttrian.bgmodel))
    tcell_to_mcell=ttrian.subfacets.facet_to_bgcell
    adapt_glue=ttrian.bgmodel.glue
    s_coarse_mfaces=adapt_glue.n2o_faces_map[end][tcell_to_mcell]
    sglue=get_glue(strian,Val(Dc))
    s_mface_to_tface=sglue.mface_to_tface
    sface_to_data_reindexed=lazy_map(Reindex(sface_to_data),s_mface_to_tface[s_coarse_mfaces])
  else 
    Gridap.Helpers.@notimplemented
  end
end 

function Gridap.Geometry.move_contributions(scell_to_val::AbstractArray,
                                            strian::GridapEmbedded.Interfaces.SubFacetTriangulation)
  model = strian.bgmodel # We cannot call get_background_model(strian) because it does not return the adapted model!!!
  ncells = num_cells(model)
  cell_to_touched = fill(false,ncells)
  scell_to_cell = strian.subfacets.facet_to_bgcell
  cell_to_touched[scell_to_cell] .= true
  Ωa = Triangulation(model,cell_to_touched)
  acell_to_val = move_contributions(scell_to_val,strian,Ωa)
  acell_to_val, Ωa 
end

function Gridap.Geometry.move_contributions(scell_to_val::AbstractArray,
                                            strian::GridapEmbedded.Interfaces.SubCellTriangulation)
  model = strian.bgmodel # We cannot call get_background_model(strian) because it does not return the adapted model!!!
  ncells = num_cells(model)
  cell_to_touched = fill(false,ncells)
  scell_to_cell = strian.subcells.cell_to_bgcell
  cell_to_touched[scell_to_cell] .= true
  Ωa = Triangulation(model,cell_to_touched)
  acell_to_val = move_contributions(scell_to_val,strian,Ωa)
  acell_to_val, Ωa 
end


# We first cut the boundary of the domain using the coarse-grid mesh 
order, ncells = 2, 8
geo = disk(0.35, x0=Point(0.5,0.5))
coarse_bgmodel = CartesianDiscreteModel((0, 1, 0, 1), (ncells, ncells))
coarse_cutgeo  = cut(coarse_bgmodel, geo)

# We then cut the boundary of the domain using the fine-grid mesh 
fine_bgmodel = Gridap.Adaptivity.refine(coarse_bgmodel, order)
fine_cutgeo  = cut(fine_bgmodel, geo)

# We modify fine_cutgeo such that all children of coarse cut cells become cut cells
# We need this to be satisfied so that we have a TestFESpace with the same number of 
# DoFs as the TrialFESpace
c2f = fine_bgmodel.glue.o2n_faces_map
c2_inoutcut = coarse_cutgeo.ls_to_bgcell_to_inoutcut[1]
for coarse_cell in 1:length(c2_inoutcut)
  if c2_inoutcut[coarse_cell] == GridapEmbedded.Interfaces.CUT
    for fine_cell in c2f[coarse_cell]
      # Inner cells stay inner cells, even if the parent is a cut cell
      # Outer cells become cut cells if the parent is a cut cell
      # if (fine_cutgeo.ls_to_bgcell_to_inoutcut[1][fine_cell]==GridapEmbedded.Interfaces.OUT)
         fine_cutgeo.ls_to_bgcell_to_inoutcut[1][fine_cell] = GridapEmbedded.Interfaces.CUT
      # end
    end
  end
end

Ω = Triangulation(fine_cutgeo, PHYSICAL)
degree = 2 * order + 1
dΩ, dΩ⁺ = Measure(Ω, degree), Measure(Ω, 12)

Ωact = Triangulation(fine_cutgeo, ACTIVE)
ΩHact = Triangulation(coarse_cutgeo, ACTIVE)

Γ = EmbeddedBoundary(fine_cutgeo)
n_Γ = get_normal_vector(Γ)
dΓ = Measure(Γ, degree)

Γg = GhostSkeleton(fine_cutgeo)
n_Γg = get_normal_vector(Γg)
dΓg = Measure(Γg, degree)

h = 1.0 / ncells  # Mesh size
γd = 10.0  # Nitsche coefficient
γg = 0.1   # ghost penalty

strategy = AggregateAllCutCells()
coarse_aggregates = aggregate(strategy,coarse_cutgeo)
fine_aggregates   = aggregate(strategy,fine_cutgeo)

colors = color_aggregates(coarse_aggregates,coarse_bgmodel)
writevtk(Triangulation(coarse_bgmodel),"trian",celldata=["aggregate"=>coarse_aggregates,"color"=>colors])

colors = color_aggregates(fine_aggregates,fine_bgmodel)
writevtk(Triangulation(fine_bgmodel),"fine_trian",celldata=["aggregate"=>fine_aggregates,"color"=>colors])


Vstd = TestFESpace(Ωact, ReferenceFE(lagrangian, Float64, 1))
Vagg = AgFEMSpace(Vstd,fine_aggregates)

Ustd = TestFESpace(ΩHact, ReferenceFE(lagrangian, Float64, order))
Uagg = TrialFESpace(AgFEMSpace(Ustd,coarse_aggregates))

#x=get_cell_dof_ids(U,Ωact)
#y=get_cell_dof_ids(U,ΩHact)
#z=get_cell_dof_ids(U,Ω)V
Gridap.get_cell_dof_ids(V,Γg)
Gridap.get_cell_is_dirichlet(U,Γg)

# V_dof_ids=get_cell_dof_ids(V,Ω)
dv=get_fe_basis(V)
du=get_trial_fe_basis(U)

∇u_Ω = Gridap.CellData.change_domain(∇(du),ReferenceDomain(),Ωact,ReferenceDomain())
∇u_Ω_n_Γg=∇u_Ω⋅n_Γg
xΓg=get_cell_points(dΓg.quad)
∇u_Ω_n_Γg.minus(xΓg)[1]

(jump(∇u_Ω⋅n_Γg)*jump(∇(dv)⋅n_Γg))(xΓg)[1]

u_Ω       = Gridap.CellData.change_domain(du,ReferenceDomain(),Ωact,ReferenceDomain())
u_Ω_plus  = u_Ω.plus(xΓg)[1][1,1]
u_Ω_minus = u_Ω..minus(xΓg)[1][1,2]


function a(u, v)
  # This is just a hack that should be replaced by a proper implementation in Gridap.jl  
  ∇u_Ω = Gridap.CellData.change_domain(∇(u),ReferenceDomain(),Ωact,ReferenceDomain())
  u_Ω  = Gridap.CellData.change_domain(u,ReferenceDomain(),Ωact,ReferenceDomain())
  #∫(∇(v)⋅∇(u)) * dΩ #+
  #∫((γd / h)*v*u - v*(n_Γ⋅∇u_Ω)-(n_Γ⋅∇(v))*u_Ω)dΓ + 
  ∫((γg * h) * (jump(u_Ω)*jump(v)))dΓg
  # ∫((γg * h) * (jump(∇u_Ω⋅n_Γg)*jump(∇(v)⋅n_Γg)))dΓg
end

l(v) =
  ∫(v * f)dΩ +
  ∫((γd / h) * v * u_ex - (n_Γ ⋅ ∇(v)) * u_ex)dΓ

dc=a(du,dv)

du=get_trial_fe_basis(U)
dv=get_fe_basis(V)
dca=a(du, dv)
dcl=l(dv)
ccm=Gridap.FESpaces.collect_cell_matrix(U, V, dca)
# ccm=Gridap.FESpaces.collect_cell_vector(V, dcl)
# ccm=Gridap.FESpaces.collect_cell_matrix_and_vector(U, V, dca, dcl)

ass=SparseMatrixAssembler(U,V)

m1 = Gridap.FESpaces.nz_counter(Gridap.FESpaces.get_matrix_builder(ass),(Gridap.FESpaces.get_rows(ass),Gridap.FESpaces.get_cols(ass)))
Gridap.FESpaces.symbolic_loop_matrix!(m1,ass,ccm)
m2 = Gridap.FESpaces.nz_allocation(m1)
Gridap.FESpaces.numeric_loop_matrix!(m2,ass,ccm)
m3 = create_from_nz(m2)
m3


A=assemble_matrix(a, U, V)
# assemble_vector(l, V)
# dc=a(du, dv)
# l(dv)

op = AffineFEOperator(a, l, U, V)
ûh = solve(op)
ûh = Gridap.CellData.change_domain(ûh,ReferenceDomain(),Ωact,ReferenceDomain())
eh = ûh - u_ex
sum(∫(eh * eh)dΩ)

# dv
# du
# dv_du=dv*du
# x=get_cell_points(Ωact)
# dv_du(x)[2][4,:,:]
# Gridap.CellData.is_change_possible(Ωact,Ω)
# get_background_model(Ω)
# This is NOT working. However, I guess it is not going to be needed, as we first transform from ΩHact to Ωact
# du_Ω = Gridap.CellData.change_domain(du, ReferenceDomain(), Ω, ReferenceDomain())
# This is working. It seems to be doing the right thing.
# dv_du_Ω = Gridap.CellData.change_domain(dv_du, ReferenceDomain(), Ω, ReferenceDomain())
# Λ  = Skeleton(fine_bgmodel)
# nΛ = get_normal_vector(Λ)
# dΛ = Measure(Λ, degree)
# t=Triangulation(fine_bgmodel)
# Vt = TestFESpace(t, ReferenceFE(lagrangian, Float64, order))
# Ut = TrialFESpace(TestFESpace(t, ReferenceFE(lagrangian, Float64, order)))
# dvt = get_fe_basis(Vt)
# dut = get_trial_fe_basis(Ut)