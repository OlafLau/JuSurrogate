using PhaseDiagram
using Documenter

makedocs(;
    modules=[PhaseDiagram],
    authors="Yi-Xin Liu <lyx@fudan.edu.cn> and contributors",
    repo="https://github.com/liuyxpp/PhaseDiagram.jl/blob/{commit}{path}#L{line}",
    sitename="PhaseDiagram.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://liuyxpp.github.io/PhaseDiagram.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/liuyxpp/PhaseDiagram.jl",
)
