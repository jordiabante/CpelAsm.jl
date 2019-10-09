push!(LOAD_PATH,"../src/")
using Documenter,CpelAsm

makedocs(sitename = "CpelAsm Documentation",
         format = Documenter.HTML(prettyurls=get(ENV,"CI",nothing)=="true"),
         authors = "Jordi Abante",
)

deploydocs(
    repo = "github.com/jordiabante/CpelAsm.jl.git",
    target = "build",
)
