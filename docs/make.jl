push!(LOAD_PATH,"../src/")
using Documenter,CpelAsm

makedocs(format = Documenter.HTML(prettyurls=get(ENV,"CI",nothing)=="true"),
         sitename = "CpelAsm.jl",
         doctest = false,
         strict = false,
         pages = [
            "Home"          => "index.md",
            "Toy Example"   => "toy_example.md"
         ],
         authors = "Jordi Abante"
)

deploydocs(
    repo = "github.com/jordiabante/CpelAsm.jl.git",
    deps = nothing,
    make = nothing,
    branch = "gh-pages"
)
