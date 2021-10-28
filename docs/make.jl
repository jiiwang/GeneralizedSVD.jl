using Documenter, GeneralizedSVD

makedocs(
    authors = "Ji Wang, Zhaojun Bai",
    sitename = "GeneralizedSVD.jl",
    modules = [GeneralizedSVD],
    pages = [
        "Home" => "index.md",
        ],
)

deploydocs(
    repo = "github.com/jiiwang/GeneralizedSVD.jl.git",
)
