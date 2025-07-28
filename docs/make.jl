using Documenter, GridapEmbedded

pages = [
    "Home" => "index.md",
    "Constructive Solid Geometry (CSG)" => "CSG.md",
    "Embedded Interfaces" => "Interfaces.md",
    "Level Set Cutters" => "LevelSetCutters.md",
    "Aggregated FEM" => "AggregatedFEM.md",
    "Moment-Fitted Quadratures" => "MomentFittedQuadratures.md",
    "Geometrical Derivatives" => "GeometricalDerivatives.md",
    "Distributed computing" => "Distributed.md",
]

makedocs(;
    modules = [GridapEmbedded],
    format = Documenter.HTML(
        size_threshold=nothing
    ),
    sitename = "GridapEmbedded.jl",
    authors = "Francesc Verdugo <f.verdugo.rojano@vu.nl>, Eric Neiva <eric.neiva@college-de-france.fr> and Santiago Badia <santiago.badia@monash.edu>",
    pages = pages,
    warnonly = false,
)

deploydocs(;
    repo="github.com/gridap/GridapEmbedded.jl",
)
