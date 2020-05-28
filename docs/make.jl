using Documenter
using LatticeRules

makedocs(;
    modules=[LatticeRules],
    authors="Pieterjan Robbe <pieterjan.robbe@kuleuven.be> and contributors",
    repo="https://github.com/PieterjanRobbe/LatticeRules.jl/blob/{commit}{path}#L{line}",
    sitename="LatticeRules.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://PieterjanRobbe.github.io/LatticeRules.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/PieterjanRobbe/LatticeRules.jl",
)
