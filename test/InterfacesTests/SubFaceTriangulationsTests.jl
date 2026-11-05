module SubFaceTriangulationsTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.CellData
using LinearAlgebra
using GridapEmbedded
using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: SubFaceTriangulation, SubCellTriangulation, SubFacetTriangulation

# Df = Dp-2: one segment inside every cell of a TET mesh, glued to the cells (Dr = Dp = 3)

bgmodel = simplexify(CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2)))
cell_map = get_cell_map(bgmodel)
ξa = Point(0.1,0.2,0.3)
ξb = Point(0.5,0.1,0.15)
ncells = num_cells(bgmodel)

cell_to_points = Table(collect(Int32,1:2ncells),collect(Int32,1:2:2ncells+1))
cell_to_bgcell = collect(Int32,1:ncells)
point_to_rcoords = repeat([ξa,ξb],ncells)
point_to_coords = [ evaluate(cell_map[cell],ξ) for cell in 1:ncells for ξ in (ξa,ξb) ]
st = SubFaceData{1}(cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords)

@test st isa SubFaceData{1,3,3}
@test !(st isa SubCellData) && !(st isa SubFacetData)
@test num_cell_dims(st) == 1
@test isnothing(st.cell_to_normal)
@test_throws Exception SubFaceData{1}(cell_to_points,cell_to_bgcell,point_to_coords,point_to_rcoords,
  [VectorValue(1.0,0.0,0.0) for i in 1:ncells])

trian = SubFaceTriangulation(st,bgmodel)
@test trian isa Triangulation{1,3}
@test !(trian isa SubCellTriangulation) && !(trian isa SubFacetTriangulation)
@test get_background_model(trian) === bgmodel
@test num_cells(trian) == ncells
@test !isnothing(get_glue(trian,Val(3)))
@test isnothing(get_glue(trian,Val(2)))
@test isnothing(get_glue(trian,Val(1)))
@test_throws Exception get_facet_normal(trian)

# A quadratic background field, integrated on the segments through the glue
f(x) = x[1] + 2*x[2]^2 - x[1]*x[3]
V = FESpace(bgmodel,ReferenceFE(lagrangian,Float64,2))
fh = interpolate(f,V)
dΛ = Measure(trian,4)
I = sum(∫(fh)dΛ)
I_exact = sum(1:ncells) do cell
  a, b = point_to_coords[2cell-1], point_to_coords[2cell]
  norm(b-a)*(f(a) + 4f((a+b)/2) + f(b))/6   # Simpson: exact for quadratics
end
@test I ≈ I_exact rtol=1e-12

# Restriction, empty, append!
st2 = SubFaceData(st,[2,5])
@test st2.cell_to_bgcell == Int32[2,5]
@test st2.cell_to_points == Table([[3,4],[9,10]])
st3 = empty(st)
@test typeof(st3) == typeof(st) && length(st3.cell_to_bgcell) == 0
append!(st3,st2); append!(st3,st2)
@test st3.cell_to_bgcell == Int32[2,5,2,5]
@test st3.cell_to_points == Table([[3,4],[9,10],[3+2ncells,4+2ncells],[9+2ncells,10+2ncells]])

# Facets without stored normals (Df = Dp-1, Dr = Dp), from the level-set cutter's subfacets

for (domain,partition,geo) in (
  ((0,1,0,1),(10,10),disk(0.31,x0=Point(0.48,0.53))),
  ((0,1,0,1,0,1),(6,6,6),sphere(0.31,x0=Point(0.48,0.53,0.51))) )

  model = CartesianDiscreteModel(domain,partition)
  D = num_cell_dims(model)
  cutgeo = cut(model,geo)

  # Aliases and former field names
  sc, sf = cutgeo.subcells, cutgeo.subfacets
  @test sc isa SubCellData{D,D} && !(sc isa SubFacetData)
  @test sf isa SubFacetData{D} && !(sf isa SubCellData)
  @test sf isa SubFaceData{D-1,D,D}
  @test sf.facet_to_normal === sf.cell_to_normal
  @test sf.facet_to_points === sf.cell_to_points
  @test sf.facet_to_bgcell === sf.cell_to_bgcell
  Ω = Triangulation(cutgeo,CUT_IN)
  Γ = EmbeddedBoundary(cutgeo)
  @test Ω isa SubCellTriangulation{D,D}
  @test Γ isa SubFacetTriangulation{D-1,D}
  @test Γ.subfacets === Γ.subfaces

  # Normals: only the stored ones, as the facet normals of the grid
  n_stored = sf.cell_to_normal
  n = length(n_stored)
  Γs = SubFaceTriangulation(sf,model)
  nΓ = get_facet_normal(Γs)
  x0 = Point(ntuple(_->0.0,D-1))
  @test all(i -> evaluate(nΓ[i],x0) == n_stored[i], 1:n)
  orientation = [ isodd(i) ? -1 : 1 for i in 1:n ]
  sf2 = SubFacetData(sf,1:n,orientation)
  @test sf2.cell_to_normal == n_stored .* orientation
  @test sf2.cell_to_points == sf.cell_to_points
  @test all(i -> evaluate(get_facet_normal(SubFaceTriangulation(sf2,model))[i],x0) == n_stored[i]*orientation[i], 1:n)

  # Facet pieces without stored normals: no normals at all
  nf = SubFaceData{D-1}(sf.cell_to_points,sf.cell_to_bgcell,sf.point_to_coords,sf.point_to_rcoords)
  @test !(nf isa SubFacetData) && !(nf isa SubCellData)
  Γn = SubFaceTriangulation(nf,model)
  @test !(Γn isa SubFacetTriangulation)
  @test_throws Exception get_facet_normal(Γn)
  @test_throws Exception SubFaceData(nf,1:n,orientation)
  @test sum(∫(1)Measure(Γn,1)) == sum(∫(1)Measure(Γs,1))

  # Active model (gluing the pieces)
  @test num_cells(get_active_model(Γn)) == n
end

# Every combination 0 ≤ Df ≤ Dr ≤ Dp (Dr ≥ 1), Dp = 2,3: one Df-simplex inside some background
# Dr-faces; a linear field integrated in physical space and through the glue (field on the
# reference Dr-faces, so the glue's reference map is used) against the exact value.

function ref_piece(Df,Dr)
  ξ0 = Point(ntuple(i->0.15,Dr))
  [ ξ0 + (k==0 ? zero(ξ0) : 0.3*Point(ntuple(i->Float64(i==k),Dr))) for k in 0:Df ]
end

flin(x) = 1.0 + x[1] + 2x[2] + (length(x)>2 ? 3x[3] : 0.0)

for Dp in 2:3, Dr in 1:Dp, Df in 0:Dr
  model = simplexify(CartesianDiscreteModel(Tuple(repeat([0,1],Dp)),Tuple(fill(2,Dp))))
  Λ = Dr == Dp ? Triangulation(model) : Triangulation(ReferenceFE{Dr},model)
  fmap = get_cell_map(Λ)
  bgfaces = collect(1:3:num_cells(Λ))
  ξs = ref_piece(Df,Dr)
  n = length(bgfaces)
  points = Table(collect(Int32,1:n*(Df+1)),collect(Int32,1:Df+1:n*(Df+1)+1))
  coords = [ evaluate(fmap[face],ξ) for face in bgfaces for ξ in ξs ]
  st = SubFaceData{Df}(points,collect(Int32,bgfaces),coords,repeat(ξs,n))
  Σ = SubFaceTriangulation(st,model)
  @test num_cell_dims(Σ) == Df && num_point_dims(Σ) == Dp
  @test isnothing(get_glue(Σ,Val(Dr == 1 ? 2 : Dr-1)))
  exact = sum(1:n) do p
    X = coords[(p-1)*(Df+1)+1:p*(Df+1)]
    G = [ (X[i+1]-X[1])⋅(X[j+1]-X[1]) for i in 1:Df, j in 1:Df ]
    vol = Df == 0 ? 1.0 : sqrt(det(G))/factorial(Df)
    vol*flin(sum(X)/(Df+1))
  end
  dΣ = Measure(Σ,1)
  fΛ = change_domain(CellField(flin,Λ),ReferenceDomain())
  @test sum(∫(flin)dΣ) ≈ exact rtol=1e-12
  @test sum(∫(fΛ)dΣ) ≈ exact rtol=1e-12
end

end # module
