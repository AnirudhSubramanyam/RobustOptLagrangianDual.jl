using Documenter, RobustOptLagrangianDual

makedocs(
    sitename="RobustOptLagrangianDual.jl",
    format=Documenter.HTML(
        prettyurls=Base.get(ENV, "CI", nothing) == "true",
        mathengine=Documenter.KaTeX()
    ),
    pages=[
        "Home" => "index.md",
        "Library" => [
            "Algorithms" => "algorithms.md",
            "Solvers" => "solvers.md",
            "Internals" => "internals.md",
            "Examples" => "examples.md",
        ],
    ]
)

deploydocs(
    devbranch="main",
    repo="github.com/AnirudhSubramanyam/RobustOptLagrangianDual.jl",
)

