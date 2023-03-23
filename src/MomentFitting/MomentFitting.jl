module MomentFitting

using Gridap
using Gridap.Arrays
using Gridap.CellData
using Gridap.Fields
using Gridap.Geometry
using Gridap.Helpers
using Gridap.Polynomials
using Gridap.ReferenceFEs

using GridapEmbedded.Interfaces
using GridapEmbedded.CSG

using FillArrays

using Gridap.Geometry: FaceToCellGlue
using Gridap.Geometry: push_normal
using Gridap.Polynomials: jacobi
using Gridap.Polynomials: _q_filter
using Gridap.Polynomials: _define_terms
using Gridap.Polynomials: _maximum
using Gridap.Polynomials: _set_value!

using GridapEmbedded.Interfaces: SubFacetTriangulation
using GridapEmbedded.Interfaces: SubFacetBoundaryTriangulation
using GridapEmbedded.Interfaces: CutInOrOut
using GridapEmbedded.Interfaces: AbstractEmbeddedDiscretization
using GridapEmbedded.Interfaces: get_geometry
using GridapEmbedded.LevelSetCutters: _check_and_get_polytope

import Gridap.Arrays: return_type
import Gridap.Arrays: return_cache
import Gridap.Arrays: evaluate!
import Gridap.Polynomials: get_order
import Gridap.Polynomials: get_orders
import Gridap.Polynomials: get_exponents

export MomentFittingMeasures
export MomentFittingQuad

include("JacobiPolynomialBases.jl")

include("CutCellMoments.jl")

end #module

