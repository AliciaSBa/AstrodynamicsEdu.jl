# Packages required for testing
using Test
using LinearAlgebra

# Package under test
using AstrodynamicsEdu

@time @testset "AstrodynamicsEdu Package Tests" begin
    testdir = joinpath(dirname(@__DIR__), "test")
    # Write your tests here.
    @time @testset "Astrodynamics.LinearAlgebraTypes" begin
        include(joinpath(testdir, "test_LinearAlgebraTypes.jl"))
        include(joinpath(testdir, "test_IdealTwoBodyProblem.jl"))
        include(joinpath(testdir, "test_OrbitalManeuvers.jl"))
        include(joinpath(testdir, "test_OrbitalPerturbations.jl"))
    end
end
