# Packages required for testing
using Test
using LinearAlgebra

# Package under test
using Astrodynamics

@time @testset "Astrodynamics Package Tests" begin
    testdir = joinpath(dirname(@__DIR__), "test")
    # Write your tests here.
    @time @testset "Astrodynamics.LinearAlgebraTypes" begin
        include(joinpath(testdir, "test_LinearAlgebraTypes.jl"))
    end
end
