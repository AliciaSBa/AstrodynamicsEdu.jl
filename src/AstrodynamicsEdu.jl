module AstrodynamicsEdu

# USINGS
using LinearAlgebra
using Plots
using Unitful
using DifferentialEquations
using Roots
using Reexport
using DelimitedFiles
using RemoteFiles

# INClUDE FILES
include("linearAlgebraTypes.jl")

include("idealTwoBodyProblem.jl")

include("orbitalPlots.jl")

include("astroConstants.jl")

include("orbitalPerturbations.jl")

include("orbitalManeuvers.jl")

end
