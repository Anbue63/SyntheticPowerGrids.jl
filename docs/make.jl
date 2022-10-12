using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(; path=dirname(@__DIR__))) # adds the package this script is called from
push!(LOAD_PATH, Base.Filesystem.abspath("../"))
using Documenter, SyntheticPowerGrids

makedocs(;
    modules=[SyntheticPowerGrids],
    authors = "Anna Büttner and contributors",
    sitename = "SyntheticPowerGrids.jl",
    pages = [
        "General" => "index.md",
        "Examples" => ["Getting started" => "easy_example.md",
                       "Mixed Network"   => "mixed_network.md",
                       "Adding new nodal dynamics" => "own_nodal_dynamics.md"]]
)