using Documenter, Cuba

makedocs(
    modules = [Cuba],
    sitename = "Cuba",
    strict = true,
)

deploydocs(
    repo = "github.com/giordano/Cuba.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
