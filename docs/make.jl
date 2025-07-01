using Documenter
using ICAforECGrecordings

makedocs(
    sitename    = "ICAforECGrecordings.jl",
    modules     = [ICAforECGrecordings],
    format      = Documenter.HTML(),
    pages       = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Usage" => "usage.md",
        "API" => "api.md",
    ],
    repo        = "https://github.com/Tim-Mueller-Bagehl/ICAforECGrecordings",
    repo_branch = "main",
)

deploydocs(github_action = true)