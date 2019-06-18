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
#=deploydocs(
    repo = "<repository url>"
)=#
