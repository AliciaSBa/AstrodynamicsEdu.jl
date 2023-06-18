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

# INCLUDE SUBMODULES
include("submodules/LinearAlgebraTypes.jl")
@reexport using .LinearAlgebraTypes

include("submodules/IdealTwoBodyProblem.jl")
@reexport using .IdealTwoBodyProblem

include("submodules/OrbitalPlots.jl")
@reexport using .OrbitalPlots

include("submodules/AstroConstants.jl")
@reexport using .AstroConstants

#include("submodules/OrbitalPerturbations.jl")
#@reexport using .OrbitalPerturbations

include("submodules/OrbitalManeuvers.jl")
@reexport using .OrbitalManeuvers

end
