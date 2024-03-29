# make.jl file
using Documenter
push!(LOAD_PATH, "C:/Users/alici/.julia/dev/AstrodynamicsEdu/src")
using AstrodynamicsEdu

#makedocs(sitename="AstrodynamicsEdu Documentation")


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
            "modules/idealTwoBodyProblem.md",
            "modules/astroConstants.md",
            "modules/orbitalManeuvers.md",
            "modules/orbitalPerturbations.md",
            "modules/orbitalPlots.md"
        ],
        "Tutorials" => Any[
            "tutorials/example1.md",
            "tutorials/example2.md",
            "tutorials/example3.md"
        ],
        # Define your documentation pages here
    ]
)

deploydocs(
    repo = "github.com/AliciaSBa/AstrodynamicsEdu.jl",
    devbranch = "main",
    devurl = "latest",
)



