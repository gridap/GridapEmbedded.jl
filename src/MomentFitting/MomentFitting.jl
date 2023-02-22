module MomentFitting

using Gridap
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields

using Gridap.CellData

using Gridap.Geometry
using Gridap.Geometry: get_reffes
using Gridap.Geometry: get_cell_to_parent_cell
using Gridap.Geometry: FaceToCellGlue
using Gridap.Geometry: push_normal

using Gridap.Polynomials
using Gridap.Polynomials: jacobi
using Gridap.Polynomials: _q_filter
using Gridap.Polynomials: _define_terms
using Gridap.Polynomials: _maximum
using Gridap.Polynomials: _set_value!

using Gridap.ReferenceFEs: get_polytope
using Gridap.ReferenceFEs: compute_nodes
using Gridap.ReferenceFEs: LagrangianDofBasis

using Gridap.Arrays: CompressedArray

using Gridap.ReferenceFEs
using Gridap.ReferenceFEs: Quadrature
using Gridap.ReferenceFEs: GenericQuadrature

using GridapEmbedded.Interfaces
using GridapEmbedded.Interfaces: SubFacetTriangulation
using GridapEmbedded.Interfaces: SubFacetBoundaryTriangulation
using GridapEmbedded.Interfaces: CutInOrOut

using GridapEmbedded.LevelSetCutters: _check_and_get_polytope

using GridapEmbedded.CSG

import Gridap.Polynomials: get_exponents
import Gridap.Polynomials: get_order
import Gridap.Polynomials: get_orders

import Gridap.Arrays: return_type
import Gridap.Arrays: return_cache
import Gridap.Arrays: evaluate!

import GridapEmbedded.Interfaces: cut_facets

import LinearAlgebra: pinv, cond
using FillArrays

export MomentFittingMeasures
export MomentFittingQuad

include("JacobiPolynomialBases.jl")

include("CutCellMoments.jl")

end #module

