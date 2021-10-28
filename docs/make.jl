# dot is needed before GSVD module, which means current Module
using Documenter, .GSVD

makedocs(
    doctest = true,
    authors = "Ji Wang, Zhaojun Bai",
    sitename = "GSVD.jl",
)
