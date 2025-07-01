using ICAforECGrecordings
using Documenter

DocMeta.setdocmeta!(ICAforECGrecordings, :DocTestSetup, :(using ICAforECGrecordings); recursive=true)

makedocs(;
    modules=[ICAforECGrecordings],
    authors="Tim Mueller-Bagehl <tim.mueller-bagehl@campus.tu-berlin.de>",
    sitename="ICAforECGrecordings.jl",
    format=Documenter.HTML(;
        canonical="https://Tim-Mueller-Bagehl.github.io/ICAforECGrecordings.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Tim-Mueller-Bagehl/ICAforECGrecordings.jl",
    devbranch="gh-pages",
)
