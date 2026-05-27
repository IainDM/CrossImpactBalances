using Documenter
using CrossImpactBalances

DocMeta.setdocmeta!(CrossImpactBalances, :DocTestSetup,
                    :(using CrossImpactBalances); recursive=true)

makedocs(
    modules = [CrossImpactBalances],
    sitename = "CrossImpactBalances.jl",
    authors = "Iain Morrow",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    pages = [
        "Home" => "index.md",
        "API reference" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/IainDM/CrossImpactBalances.jl",
    devbranch = "main",
    push_preview = true,
)
