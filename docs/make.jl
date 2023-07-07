# make.jl file
using Documenter
push!(LOAD_PATH,"C:/Users/alici/.julia/dev/AstrodynamicsEdu/src")
using AstrodynamicsEdu

makedocs(
    name = "YourPackageName",
    pages = [
        "Home" => "index.md",
        # Add more pages if needed
    ],
    format = Documenter.HTML()
)

makedocs(
    modules   = [AstrodynamicsEdu],  # Replace with the modules from your package
    doctest   = false,
    clean     = true,
    linkcheck = true,
    logo = "assets/logo.png",
    format    = Documenter.HTML(),
    sitename  = "AstrodynamicsEdu.jl",
    authors   = "Alicia Sanjurjo",
    pages     = Any[
        "Home" => "index.md",
        "Modules" => Any[
            "modules/linearAlgebraTypes.md",
            "modules/idealTwoBodyProblem.md"
        ],
        "Tutorials" => Any[
            "tutorials/example1.md",
        ],
        # Define your documentation pages here
    ]
)

deploydocs(
    repo = "github.com/AliciaSBa/AstrodynamicsEdu.jl",
    devbranch = "main",
    devurl = "latest",
)

