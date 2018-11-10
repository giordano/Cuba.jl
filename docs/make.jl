using Documenter, Cuba

makedocs(
    modules = [Cuba],
    format = :html,
    sitename = "Cuba",
)

deploydocs(
    repo = "github.com/giordano/Cuba.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
