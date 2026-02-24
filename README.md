# GridapEmbedded <img src="https://github.com/gridap/Gridap.jl/blob/master/images/color-logo-only.png" width="40" title="Gridap logo">

Embedded finite element methods, level set surface descriptions and constructive solid geometry.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/GridapEmbedded.jl/stable)
[![Build Status](https://github.com/gridap/GridapEmbedded.jl/workflows/CI/badge.svg?branch=master)](https://github.com/gridap/GridapEmbedded.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/gridap/GridapEmbedded.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/GridapEmbedded.jl)

## Installation

```julia
# Type ] to enter package mode
pkg> add GridapEmbedded 
```

## Examples

### Constructive Solid Geometry (CSG)

```julia
julia> include("examples/PoissonCSGCutFEM/PoissonCSGCutFEM.jl")
julia> PoissonCSGCutFEM.main(n=40,outputfile="results1")
```

<img src="https://upload.wikimedia.org/wikipedia/commons/8/8b/Csg_tree.png" width="300"><img src="https://github.com/gridap/GridapEmbedded.jl/blob/master/examples/PoissonCSGCutFEM/PoissonCSGCutFEM_solution.png?raw=true" width="300">

*left picture by wikipedia.org*

```julia
julia> include("examples/StokesTubeWithObstacleCutFEM/StokesTubeWithObstacleCutFEM.jl")
julia> StokesTubeWithObstacleCutFEM.main(n=10,outputfile="results2")
```

<img src="https://github.com/gridap/GridapEmbedded.jl/blob/master/examples/StokesTubeWithObstacleCutFEM/StokesTubeWithObstacleCutFEM_solution.png?raw=true" width="600">

### Bimaterial problems

```julia
julia> include("examples/BimaterialLinElastCutFEM/BimaterialLinElastCutFEM.jl")
julia> BimaterialLinElastCutFEM.main(n=4,outputfile="results3")
```

<img src="https://raw.githubusercontent.com/gridap/GridapEmbedded.jl/master/examples/BimaterialLinElastCutFEM/BimaterialLinElastCutFEM_solution.png" width="400">
