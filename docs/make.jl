using MotivicHomotopy
using Documenter

DocMeta.setdocmeta!(MotivicHomotopy, :DocTestSetup, :(using MotivicHomotopy); recursive=true)

makedocs(;
    modules=[MotivicHomotopy],
    authors="Stephanie Atherton, Nikita Borisov, Thomas Brazelton, Somak Dutta, Frenly Espino, Thomas Hagedorn, Zhaobo Han, Jordy Lopez Garcia, Joel Louwsma, Yuyuan Luo, Wern Juin Gabriel Ong, Ruzho Sagayaraj, Andrew Tawfeek",
    sitename="MotivicHomotopy.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaMotivicHomotopy.github.io/MotivicHomotopy.jl",
        edit_link="main",
        sidebar_sitename=false,
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API reference" => "api.md",
    ],
    # The examples in docstrings are captured REPL sessions, not run at build
    # time (they require Oscar to be loaded and exact display matching).
    doctest=false,
)

deploydocs(;
    repo="github.com/JuliaMotivicHomotopy/MotivicHomotopy.jl",
    devbranch="main",
)
