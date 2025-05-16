
"""
    abstract type Cutter <: GridapType end

Abstract type for all mesh cutters. Has to be paired with a [`CSG.Geometry`](@ref) to 
cut the background mesh.

## Methods

- [`cut(cutter::Cutter,background,geom)`](@ref)
- [`cut_facets(cutter::Cutter,background,geom)`](@ref)
- [`compute_bgcell_to_inoutcut(cutter::Cutter,background,geom)`](@ref)
- [`compute_bgfacet_to_inoutcut(cutter::Cutter,background,geom)`](@ref)

Generally `cut` and `cut_facets` dispatch based on the geometry provided, so it is 
generally more convennient to call the following methods instead:

- [`cut(background,geom)`](@ref)
- [`cut_facets(background,geom)`](@ref)

"""
abstract type Cutter <: GridapType end

"""
    cut(cutter::Cutter,background,geom)

Cut the background mesh with the provided cutter and geometry, returnning the cut cells.
The cut cells are returned as an [`EmbeddedDiscretization`](@ref) object.
"""
function cut(cutter::Cutter,background,geom)
  @abstractmethod
end

"""
    compute_bgcell_to_inoutcut(cutter::Cutter,background,geom)

Returns an array of IN/OUT/CUT states for each cell in the background mesh.
"""
function compute_bgcell_to_inoutcut(cutter::Cutter,background,geom)
  @abstractmethod
end

"""
    cut_facets(cutter::Cutter,background,geom)

Cut the background mesh with the provided cutter and geometry, returning the cut facets.
The cut facets are returned as an [`EmbeddedFacetDiscretization`](@ref) object.
"""
function cut_facets(cutter::Cutter,background,geom)
  @abstractmethod
end

"""
    compute_bgfacet_to_inoutcut(cutter::Cutter,background,geom)

Returns an array of IN/OUT/CUT states for each facet in the background mesh.
"""
function compute_bgfacet_to_inoutcut(cutter::Cutter,background,geom)
  @abstractmethod
end

function EmbeddedDiscretization(cutter::Cutter,background,geom)
  cut(cutter,background,geom)
end

function EmbeddedFacetDiscretization(cutter::Cutter,background,geom)
  cut_facets(cutter,background,geom)
end
