using SpgpGpu
using Documenter

DocMeta.setdocmeta!(SpgpGpu, :DocTestSetup, :(using SpgpGpu); recursive=true)

makedocs(;
    modules=[SpgpGpu],
    authors="Deepak Selvan <dselvan@gmail.com> and contributors",
    repo="https://github.com/dselvan/SpgpGpu.jl/blob/{commit}{path}#{line}",
    sitename="SpgpGpu.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dselvan.github.io/SpgpGpu.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dselvan/SpgpGpu.jl",
)
