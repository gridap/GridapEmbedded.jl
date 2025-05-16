
# Embedded Interfaces

```@meta
CurrentModule = GridapEmbedded.Interfaces
```

## Domain Nomenclature

Throughout this documentation, many methods accept arguments that select different parts of the cut domain. We split the domain into the following parts:

**The background mesh entities** (cells, facets, nodes) are classified as `IN`, `OUT` or `CUT`. The `IN` and `OUT` background cells are uncut, i.e completely inside or outside the geometry, respectively. These states are internally defined as constants:

```julia
  const IN = -1
  const OUT = 1
  const CUT = 0
```

**The `CUT` background cells** are cut by the embedded boundary, and split into subcells/subfacets. The subcells/subfacets are classified as `IN` or `OUT` depending on whether they are inside or outside the geometry. `CUT_IN` and `CUT_OUT` subentities can be accessed using the `CutInOrOut` objects:

```julia
  struct CutInOrOut
    in_or_out::Int
  end
  const CUT_IN = CutInOrOut(IN)
  const CUT_OUT = CutInOrOut(OUT)
```

For FEM, we generally want to get sets of uncut and cut cells together, for a given state `IN/OUT`. These are referred as `PHYSICAL` parts of the domain. Moreover, FE spaces are generally defined over the background mesh and need to span both `IN/OUT` and `CUT` background cells. These are referred as `ACTIVE` parts of the domain. You can extract the `PHYSICAL` and `ACTIVE` parts of the domain using the following constants:

```julia
const PHYSICAL_IN = (CUT_IN,IN)
const PHYSICAL_OUT = (CUT_OUT,OUT)
const PHYSICAL = PHYSICAL_IN
```

```julia
struct ActiveInOrOut
  in_or_out::Int
end
const ACTIVE_IN = ActiveInOrOut(IN)
const ACTIVE_OUT = ActiveInOrOut(OUT)
const ACTIVE = ACTIVE_IN
```

## Cutters

Cutters are used to cut the background mesh according to a provided geometry.

```@autodocs
Modules = [Interfaces,]
Order   = [:type, :constant, :macro, :function]
Pages   = [
  "/Cutters.jl", 
]
```

We provide several types of cutters, including:

- **Level-Set Cutters**: Cutters for Level-Set and function-defined geometries. See [Level-Set Cutters](@ref).
- **STL Cutters**: Cutters for STL based geometries. Provided by [STLCutters.jl](https://github.com/gridap/STLCutters.jl).

## Embedded Discretizations

After cutting the background mesh, you will be returned an `EmbeddedDiscretization` object. These contain all the information you need to generate the integration meshes for embedded methods.

```@docs
AbstractEmbeddedDiscretization
EmbeddedDiscretization
EmbeddedFacetDiscretization
```

## Embedded Triangulations

From `EmbeddedDiscretization` objects, you can extract all the triangulations you need to perform integration for embedded methods. We currently provide the following methods:

```@docs
Gridap.Geometry.Triangulation(::EmbeddedDiscretization,::Any)
EmbeddedBoundary(::EmbeddedDiscretization)
GhostSkeleton(::EmbeddedDiscretization)
Gridap.Geometry.BoundaryTriangulation(::EmbeddedFacetDiscretization,::Any)
Gridap.Geometry.SkeletonTriangulation(::EmbeddedFacetDiscretization,::Any)
SubCellTriangulation
SubFacetTriangulation
SubFacetBoundaryTriangulation
```
