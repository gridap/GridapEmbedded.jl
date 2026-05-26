# Tests for the fix to issue #115:
#   https://github.com/gridap/GridapEmbedded.jl/issues/115
#
# When a compound CSG geometry (e.g. setdiff) is used, the conservative cell-level
# classification can mark background cells as CUT even when no physical subcell of
# the target domain lies in them. Before the fix, Triangulation(ACTIVE) and
# GhostSkeleton would include these phantom cells, producing isolated elements.

module Issue115Tests

using Test
using Gridap
using Gridap.Geometry
using Gridap.Adaptivity
using GridapEmbedded
using GridapEmbedded.LevelSetCutters

n = 12
base_model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(n,n)))
ref_model  = refine(base_model, refinement_method = "barycentric")
model      = get_model(ref_model)

# Two concentric circles: φ1 (radius r1=0.3), φ2 (radius r2 < r1 → empty setdiff)
r1 = 0.3
r2 = r1 - 0.5/n  # r2 < r1: Ω2 = setdiff(inside φ2, inside φ1) is geometrically empty
φ1((x,y)) = (x-0.5)^2 + (y-0.5)^2 - r1^2
φ2((x,y)) = (x-0.5)^2 + (y-0.5)^2 - r2^2

V_φ = TestFESpace(model, ReferenceFE(lagrangian,Float64,1))
φh1 = interpolate(φ1, V_φ)
φh2 = interpolate(φ2, V_φ)

geo1   = DiscreteGeometry(φh1, model, name="geo1")
geo2   = DiscreteGeometry(φh2, model, name="geo2")
Ω1_geo = setdiff(!geo1, geo2, name="Ω1")  # ring outside φ2 and inside φ1 (non-empty)
Ω2_geo = setdiff(geo2,  geo1, name="Ω2")  # inside φ2 and outside φ1 (empty, r2 < r1)
cutgeo = cut(model, union(Ω1_geo, Ω2_geo))

# Empty domain: ACTIVE triangulation must contain no cells
Ω2_act = Triangulation(cutgeo, ACTIVE, "Ω2")
@test num_cells(Ω2_act) == 0

# Empty domain: GhostSkeleton must contain no facets
Λ2 = GhostSkeleton(cutgeo, ACTIVE_IN, "Ω2")
@test num_cells(Λ2) == 0

# Non-empty domain Ω1 must be unaffected
Ω1_act = Triangulation(cutgeo, ACTIVE, "Ω1")
@test num_cells(Ω1_act) > 0

Λ1 = GhostSkeleton(cutgeo, ACTIVE_IN, "Ω1")
@test num_cells(Λ1) > 0

end # module
