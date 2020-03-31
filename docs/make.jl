using Documenter,CpelAsm

makedocs(
    sitename="CpelAsm",
    pages = [
        "Home"          => "index.md",
        "Main Commands"   => "main_commands.md",
        "Toy Example"   => "toy_example.md"
    ],
    authors = "Jordi Abante"
)

deploydocs(
    repo = "github.com/jordiabante/CpelAsm.jl.git",
)