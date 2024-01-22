using Documenter, GridapEmbedded

makedocs(;
    modules=[GridapEmbedded],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/gridap/GridapEmbedded.jl/blob/{commit}{path}#L{line}",
    sitename="GridapEmbedded.jl",
    authors="Francesc Verdugo <f.verdugo.rojano@vu.nl>, Eric Neiva <eric.neiva@college-de-france.fr> and Santiago Badia <santiago.badia@monash.edu>",
)

deploydocs(;
    repo="github.com/gridap/GridapEmbedded.jl",
)
