# GridapEmbedded

Embedded finite element methods, level set surface descriptions and constructive solid geometry.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/GridapEmbedded.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/GridapEmbedded.jl/dev)
[![Build Status](https://travis-ci.com/gridap/GridapEmbedded.jl.svg?branch=master)](https://travis-ci.com/gridap/GridapEmbedded.jl)
[![Codecov](https://codecov.io/gh/gridap/GridapEmbedded.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/GridapEmbedded.jl)

## Constructive Solid Geometry (CSG)

<img src="https://upload.wikimedia.org/wikipedia/commons/8/8b/Csg_tree.png" width="300"> 

*source: wikipedia.org*

```julia
julia> include("examples/PoissonCSGCutFEM/PoissonCSGCutFEM.jl")
julia> PoissonCSGCutFEM.main(n=40,outputfile="results1")
```

<img src="https://github.com/gridap/GridapEmbedded.jl/blob/preparing_release/examples/PoissonCSGCutFEM/PoissonCSGCutFEM_solution.png?raw=true" width="300">



