using Documenter, GridapEmbedded

makedocs(;
    modules=[GridapEmbedded],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/gridap/GridapEmbedded.jl/blob/{commit}{path}#L{line}",
    sitename="GridapEmbedded.jl",
    authors="Francesc Verdugo <fverdugo@cimne.upc.edu>",
)

deploydocs(;
    repo="github.com/gridap/GridapEmbedded.jl",
)
