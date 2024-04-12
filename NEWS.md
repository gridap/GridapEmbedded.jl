# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.9.1] - 2024-04-12

### Changed

- Update compats to include Gridap 0.18. Since PR [#80](https://github.com/gridap/GridapEmbedded.jl/pull/80).
- Optimise computation of normal displacement in dynamic simulations. Since PR [#78](https://github.com/gridap/GridapEmbedded.jl/pull/78).

### Fixed

- Compute signed distance FE function in 3D. Since PR [#77](https://github.com/gridap/GridapEmbedded.jl/pull/77).

## [0.9.0] - 2024-01-22

### Added
- Interface to Algoim v0.2, which provides algoim's high-order quadrature algorithms for domains implicitly-defined by multivariate polynomials and high-order accurate algorithms for computing closest points on implicitly-defined surfaces. Since PR [#76](https://github.com/gridap/GridapEmbedded.jl/pull/76).

## [0.8.3] - 2024-01-02

### Added
- Support for FillArrays v1. Since PR [#75](https://github.com/gridap/GridapEmbedded.jl/pull/75).
- Support for AbstractTrees v0.4. Since PR [#73](https://github.com/gridap/GridapEmbedded.jl/pull/73).

## [0.8.2] - 2023-08-22

### Added
- Moment fitted machinery. Since PR [#68](https://github.com/gridap/GridapEmbedded.jl/pull/68).

## [0.8.0] - 2021-11-03

### Changed

 - Switch to Gridap 0.17. This required some changes in the user API. In particular, one needs to specify now if you want to build a triangulation of the "physical" domain or of the "active" domain. E.g, `Triangulation(cutgeo,PHYSICAL)` or `Triangulation(cutgeo,ACTIVE)`. Taking advantage that one can define FE spaces on "body-fitted" triangulation in Gridap v0.17, the user builds the FE space in embedded simulation by using the "active" triangulation.

## [0.7.0] - 2021-06-29

### Changed

 - Switch to Gridap 0.16. Note that we specifically need the fixes introduced in Gridap v0.16.3.

## [0.6.1] - 2021-06-24

### Added
- Method `AnalyticalGeometry(::Function)` to simplify the setup of geometries from user-defined functions.
- `export AnalyticalGeometry` in `GridapEmbedded`.

## [0.6.0] - 2021-01-12

### Changed

 - Switch to Gridap 0.15.
 - Switch from TravisCI to GitHubCI

## [0.5.0] - 2020-10-21

### Added

 - `AggregateCutCellsByThreshold` strategy as a general case of `AggregateAllCutCells`. Basic usage:

   ```
     strategy = AggregateAllCutCells() # Same as before
     strategy = AggregateCutCellsByThreshold(0.5) # New
   ```

### Changed

 - Signature of methods that implement aggregation strategies

## [0.4.0] - 2020-09-19

### Added

 - Support for Gridap v0.14.
 - Optional kw-argument `edges` to `square` function.

## [0.3.1] - 2020-07-25

### Added

 - Support for Gridap v0.12 and v0.13.
 - `square` function.

## [0.3.0] - 2020-06-23

### Changed

 - Updated code according to Gridap v0.11.

## [0.2.0] - 2020-05-21

### Added

 - Extended the `Cutter` interface with function `cut_facets` that allows one to
 cut the facets of the background model instead of the cells.
 - Support for AgFEM methods. This includes the generation of cell aggregates
 (via the `aggregate` function) and the construction of aggregated FE spaces
 (via the constructor `AgFEMSpace`).

### Changed

 - Minor change in type signature of `Subtriangulation`.

## [0.1.0] - 2020-04-21

### Added

- CSG machinery
- The `Cutter` interface representing a method that "cuts" the
  cells of a background `DiscreteModel` in order to generate different types of
  `Triangulations` that allow one to integrate on cut cells and embedded boundaries.
- Concrete implementation `LevelSetCutter <: Cutter` that allows one to cut a background model
  with either an analytical levelset function, i.e. a `AnalyticalGeometry` object, or a discrete
  levelset, i.e., a `DiscreteGeometry` object.

