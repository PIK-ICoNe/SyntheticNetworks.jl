using Documenter
using SyntheticNetworks

makedocs(
    sitename = "SyntheticNetworks",
    format = Documenter.HTML(),
    modules = [SyntheticNetworks],
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

deploydocs(repo = "github.com/luap-pik/SyntheticNetworks.git")
