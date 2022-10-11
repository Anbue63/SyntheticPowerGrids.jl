using Documenter, SyntheticPowerGrids

makedocs(;
    modules=[SyntheticPowerGrids],
    authors = "Anna Büttner and contributors",
    sitename = "SyntheticPowerGrids.jl",
    pages = [
        "SyntheticPowerGrids Docs" => "index.md"]
)