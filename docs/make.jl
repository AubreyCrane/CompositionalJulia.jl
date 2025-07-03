using CompositionalJulia
using Documenter

DocMeta.setdocmeta!(CompositionalJulia, :DocTestSetup, :(using CompositionalJulia); recursive=true)

makedocs(;
    modules=[CompositionalJulia],
    authors="A. William Crane",
    sitename="CompositionalJulia.jl",
    format=Documenter.HTML(;
        canonical="https://AubreyCrane.github.io/CompositionalJulia.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/AubreyCrane/CompositionalJulia.jl",
    devbranch="master",
)
