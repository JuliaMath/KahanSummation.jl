using KahanSummation, Documenter

makedocs(
    modules = [KahanSummation],
    clean = false,
    format = :html,
    sitename = "KahanSummation.jl",
    authors = "Jeff Bezanson, Jeffrey Sarnoff, and other contributors",
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(
    julia = "nightly",
    repo = "github.com/JuliaMath/KahanSummation.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
