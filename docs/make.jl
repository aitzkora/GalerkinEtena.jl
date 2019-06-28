using Documenter
using GalerkinEtena

makedocs(
    sitename = "GalerkinEtena",
    format = Documenter.HTML(),
    modules = [GalerkinEtena]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    deps = Deps.pip("mkdocs", "python-markdown-math")
    repo = "github.com/aitzkora/GalerkinEtena.jl.git"
)
