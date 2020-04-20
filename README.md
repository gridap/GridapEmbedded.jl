# GridapEmbedded

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://gridap.github.io/GridapEmbedded.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://gridap.github.io/GridapEmbedded.jl/dev)
[![Build Status](https://travis-ci.com/gridap/GridapEmbedded.jl.svg?branch=master)](https://travis-ci.com/gridap/GridapEmbedded.jl)
[![Codecov](https://codecov.io/gh/gridap/GridapEmbedded.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/gridap/GridapEmbedded.jl)


## Example

Solve the Poisson equation `-div(grad(u)) = f` on a geometrical component obtained via constructive solid geometry (CSG). The problem is solved with an embedded finite element method using the tools provided by this package (no unstructured grid is needed).

We consider the pice used in [wikipedia to illustrate CSG.](https://en.wikipedia.org/wiki/Constructive_solid_geometry)

<img src="https://upload.wikimedia.org/wikipedia/commons/8/8b/Csg_tree.png" width="300"> 

*source: wikipedia.org*

We impose a dirichlet boundary condition `ud = 1` on the walls of the internal cylinders. Homogeneous Neumann conditions are imposed elsewehere on the boundary and a source of `f = 10` is imposed in the bulk of the piece. The obtained  solution is depicted as follows.

<img src="https://github.com/gridap/GridapEmbedded.jl/blob/preparing_release/examples/PoissonCSGCutFEM_solution.png?raw=true" width="300"> 

Find the source code to reproduce this example in this file: [examples/PoissonCSGCutFEM.jl](https://github.com/gridap/GridapEmbedded.jl/blob/preparing_release/examples/PoissonCSGCutFEM.jl)
