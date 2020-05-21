abstract type Cutter <: GridapType end

"""
    cut(cutter::Cutter,background,geom) -> EmbeddedDiscretization
"""
function cut(cutter::Cutter,background,geom)
  @abstractmethod
end

function compute_bgcell_to_inoutcut(cutter::Cutter,background,geom)
  @abstractmethod
end

"""
    cut_facets(cutter::Cutter,background,geom) -> EmbeddedDiscretization
"""
function cut_facets(cutter::Cutter,background,geom)
  @abstractmethod
end

function compute_bgfacet_to_inoutcut(cutter::Cutter,background,geom)
  @abstractmethod
end

function EmbeddedDiscretization(cutter::Cutter,background,geom)
  cut(cutter,background,geom)
end

function EmbeddedFacetDiscretization(cutter::Cutter,background,geom)
  cut_facets(cutter,background,geom)
end
