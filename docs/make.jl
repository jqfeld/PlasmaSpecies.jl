using PlasmaSpecies
using Documenter

DocMeta.setdocmeta!(PlasmaSpecies, :DocTestSetup, :(using PlasmaSpecies); recursive=true)

makedocs(;
    modules=[PlasmaSpecies],
    authors="Jan Kuhfeld <jan.kuhfeld@rub.de> and contributors",
    repo="https://github.com/jqfeld/PlasmaSpecies.jl/blob/{commit}{path}#{line}",
    sitename="PlasmaSpecies.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jqfeld.github.io/PlasmaSpecies.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jqfeld/PlasmaSpecies.jl",
    devbranch="main",
)
