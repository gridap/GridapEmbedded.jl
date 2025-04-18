using Documenter, GridapEmbedded

pages = [
    "Home" => "index.md",
    "Interfaces" => "Interfaces.md",
    "Constructive Solid Geometry (CSG)" => "CSGCutters.md",
    "Level Sets" => "LevelSetCutters.md",
    "Aggregated FEM" => "AggregatedFEM.md",
    "Moment-Fitted Quadratures" => "MomentFittedQuadratures.md",
    "Distributed computing" => "Distributed.md",
],

makedocs(;
    modules = [GridapEmbedded],
    format = Documenter.HTML(
        size_threshold=nothing
    ),
    sitename = "GridapEmbedded.jl",
    authors = "Francesc Verdugo <f.verdugo.rojano@vu.nl>, Eric Neiva <eric.neiva@college-de-france.fr> and Santiago Badia <santiago.badia@monash.edu>",
    pages = pages,
    warnonly = true,
)

deploydocs(;
    repo="github.com/gridap/GridapEmbedded.jl",
)
