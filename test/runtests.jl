# Packages required for testing
using Test
using LinearAlgebra
using DifferentialEquations

# Package under test
using AstrodynamicsEdu

@time @testset "AstrodynamicsEdu Package Tests" begin
    testdir = joinpath(dirname(@__DIR__), "test")
    # Write your tests here.
    @time @testset "AstrodynamicsEdu.linearAlgebraTypes" begin
        include(joinpath(testdir, "test_LinearAlgebraTypes.jl"))
    end
    @time @testset "AstrodynamicsEdu.idealTwoBodyProblem" begin
        include(joinpath(testdir, "test_IdealTwoBodyProblem.jl"))
    end
    @time @testset "AstrodynamicsEdu.orbitalManeuvers" begin
        include(joinpath(testdir, "test_OrbitalManeuvers.jl"))
    end
    @time @testset "AstrodynamicsEdu.orbitalPerturbations" begin
        include(joinpath(testdir, "test_OrbitalPerturbations.jl"))
    end
end
