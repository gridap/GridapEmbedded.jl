# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2020-6-23

### Changed

 - Updated code according to Gridap v0.11.

## [0.2.0] - 2020-5-21

### Added

 - Extended the `Cutter` interface with function `cut_facets` that allows one to
 cut the facets of the background model instead of the cells.
 - Support for AgFEM methods. This includes the generation of cell aggregates
 (via the `aggregate` function) and the construction of aggregated FE spaces
 (via the constructor `AgFEMSpace`).

### Changed

 - Minor change in type signature of `Subtriangulation`.

## [0.1.0] - 2020-4-21

### Added

- CSG machinery
- The `Cutter` interface representing a method that "cuts" the
  cells of a background `DiscreteModel` in order to generate different types of
  `Triangulations` that allow one to integrate on cut cells and embedded boundaries.
- Concrete implementation `LevelSetCutter <: Cutter` that allows one to cut a background model
  with either an analytical levelset function, i.e. a `AnalyticalGeometry` object, or a discrete
  levelset, i.e., a `DiscreteGeometry` object.

