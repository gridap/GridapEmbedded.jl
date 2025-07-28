module GridapEmbedded

include("CSG/CSG.jl")

include("Interfaces/Interfaces.jl")

include("LevelSetCutters/LevelSetCutters.jl")

include("AgFEM/AgFEM.jl")

include("MomentFittedQuadratures/MomentFittedQuadratures.jl")

include("AlgoimUtils/AlgoimUtils.jl")

include("Distributed/Distributed.jl")

include("Exports.jl")

include("BGP/BGP.jl")

end # module
