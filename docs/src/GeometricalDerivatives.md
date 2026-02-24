# Geometrical Derivatives

```@meta
CurrentModule = GridapEmbedded
```

The geometrical differentiation capabilities are based on the following work:

!!! note "Reference"
    "Level-set topology optimisation with unfitted finite elements and automatic shape differentiation",
    by Z. J. Wegert, J. Manyer, C. Mallon, S. Badia, V. J. Challis (2025)

To see examples of usage, please refer to the tests in `test/LevelSetCuttersTests/GeometricalDifferentiationTests.jl`.

## Discretize then differentiate

```@autodocs
Modules = [GridapEmbedded.Interfaces]
Order   = [:type, :constant, :macro, :function]
Pages   = [
  "/CutFaceBoundaryTriangulations.jl", 
]
```

## Autodiff

```@autodocs
Modules = [GridapEmbedded.LevelSetCutters]
Order   = [:type, :constant, :macro, :function]
Pages   = [
  "/DifferentiableTriangulations.jl", 
]
```
