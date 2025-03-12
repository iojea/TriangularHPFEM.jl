using TriangularhpFEM
using Documenter

DocMeta.setdocmeta!(TriangularhpFEM, :DocTestSetup, :(using TriangularhpFEM); recursive=true)

makedocs(;
    modules=[TriangularhpFEM],
    authors="Ignacio Ojea",
    sitename="TriangularhpFEM.jl",
    format=Documenter.HTML(;
        canonical="https://iojea.github.io/TriangularhpFEM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/iojea/TriangularhpFEM.jl",
    devbranch="main",
)
