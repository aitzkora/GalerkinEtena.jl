using Documenter
using GalerkinEtena

makedocs(
    sitename = "GalerkinEtena",
    authors = "Fuentes Marc",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
                )
              ))
         ),
    modules = [GalerkinEtena],
    pages = ["Documentation" => "index.md"]
        )

deploydocs(
    repo = "github.com/aitzkora/GalerkinEtena.jl.git"
    )
