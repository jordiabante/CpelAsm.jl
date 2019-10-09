push!(LOAD_PATH,"../src/")
using Documenter,CpelAsm

makedocs(sitename="CpelAsm Documentation")

deploydocs(
    repo = "github.com/jordiabante/CpelAsm.jl.git",
)
