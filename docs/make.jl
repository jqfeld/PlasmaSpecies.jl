using PlasmaSpecies
using Documenter

DocMeta.setdocmeta!(PlasmaSpecies, :DocTestSetup, :(using PlasmaSpecies); recursive=true)

makedocs(;
    modules=[PlasmaSpecies],
    authors="Jan Kuhfeld <jan.kuhfeld@rub.de> and contributors",
    repo=Remotes.GitHub("jqfeld", "PlasmaSpecies.jl"),
    sitename="PlasmaSpecies.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://jqfeld.github.io/PlasmaSpecies.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Species" => "species.md",
        "Reactions" => "reactions.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/jqfeld/PlasmaSpecies.jl",
    devbranch="main",
)
