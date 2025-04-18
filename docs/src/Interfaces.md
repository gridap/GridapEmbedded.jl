
# Embedded Interfaces

```@meta
CurrentModule = GridapEmbedded.Interfaces
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

- **CSG Cutters**: Constructive Solid Geometry (CSG) Cutters. See [Constructive Solid Geometry (CSG) Cutters](@ref).
- **Level-Set Cutters**: Level-Set Cutters. See [Level-Set Cutters](@ref).
- **STL Cutters**: Cutters for STL based geometries. Provided by [STLCutters.jl](https://github.com/gridap/STLCutters.jl).

## Embedded Discretizations

After cutting the background mesh, you will be returned an `EmbeddedDiscretization` object. These contain all the information you need to generate the integration meshes for embedded methods.

```@autodocs
Modules = [Interfaces,]
Order   = [:type, :constant, :macro, :function]
Pages   = [
  "/EmbeddedDiscretizations.jl", 
  "/EmbeddedFacetDiscretizations.jl", 
]
```

## Embedded Triangulations

```@autodocs
Modules = [Interfaces,]
Order   = [:type, :constant, :macro, :function]
Pages   = [
  "/SubCellTriangulations", 
  "/SubFacetTriangulations.jl"
]
```
