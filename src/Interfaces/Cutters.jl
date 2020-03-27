abstract type Cutter <: GridapType end

"""
    cut(cutter::Cutter,background,geom) -> EmbeddedDiscretization
"""
function cut(cutter::Cutter,background,geom)
  @abstractmethod
end

function EmbeddedDiscretization(cutter::Cutter,background,geom)
  cut(cutter,background,geom)
end
