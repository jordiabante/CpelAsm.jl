push!(LOAD_PATH,"../src/")
using Documenter,CpelAsm

makedocs(sitename = "CpelAsm.jl",
         format = Documenter.HTML(prettyurls=get(ENV,"CI",nothing)=="true"),
         authors = "Jordi Abante",
         pages = ["Home" => "index.md"],
)

deploydocs(
    repo = "github.com/jordiabante/CpelAsm.jl.git",
    deps = nothing,
    make = nothing
)
