# ORBITAL PERTURBATIONS TEST

#=
using LinearAlgebra
include("linearAlgebraTypes.jl")
include("idealTwoBodyProblem.jl")
include("orbitalPerturbations.jl")
inlcude("astroConstants.jl")
using DifferentialEquations
using Test
=#

# Test the functions in the OrbitalPerturbations module

# Test the J2 acceleration function
function test_J2_acceleration()
    rr = [-5200,-3500, 4637] # km/s
    ap_J2 = J2_acceleration(rr)
    J2_expected = [-3.643729e-06,-2.452510e-06,-5.207047e-06]
    @test isapprox(ap_J2, J2_expected, rtol=1e-6, atol=1e-8) # Expected acceleration
end

# Test the drag acceleration function
function test_drag_acceleration()
    rr = [-5200,-3500, 4637] # km/s
    r = norm(rr) # km
    zeta = r - 6371.0 # km
    vv = [5.2,6.4,2.5] # km/s
    bc = 2.0*10^6 # kg/km^2
    ap_drag = drag_acceleration(zeta, vv, bc)
    expected_drag = [-3.408820e-16,-4.195470e-16,-1.638856e-16] # km/s^2
    @test isapprox(ap_drag, expected_drag, rtol=1e-6, atol=1e-8) # Expected acceleration
end

# Test the cowell function
function test_cowell()
    r0 = [-5200.0,-3500.0, 4637.0] # km/s
    v0 = [5.2,6.4,2.5] # km/s
    Bc = 2.0*10^6 # kg/km^2
    t = [0.0, 2*30*24*3600.0, 10] # s (2 months)
    # 1. No perturbations case
    flag_drag = false 
    flag_J2 = false 
    r1, v1 = cowell(r0, v0, Bc, t, flag_drag, flag_J2)
    println("r1= ", r1)
    # Define expected results for the test case
    expected_r1 = [-5200 4.043618e+03 -2.414651e+03 -8.120478e+03;
                        -3500 -1.504363e+03 -1.046909e+04 -1.505092e+04;
                         4637 -1.400622e+04 -1.961179e+04 -1.634826e+04]
    expected_v1 = [5.200000e+00 -1.567788e+00 -1.903140e+00 -1.337834e+00;
                        6.400000e+00 -3.146063e+00 -2.006148e+00 -6.225777e-01; 
                        2.500000e+00 -3.747539e+00 -8.762529e-02 1.876878e+00]
    @test isapprox(r1, expected_r1, rtol=1e-6, atol=1e-8) # Compare calculated and expected position vectors
    @test isapprox(v1, expected_v1, rtol=1e-6, atol=1e-8) # Compare calculated and expected velocity vectors
    # 2. Drag perturbation case
    flag_drag = true 
    flag_J2 = false
    r2, v2 = cowell(r0, v0, Bc, t, flag_drag, flag_J2)
    # Define expected results for the test case
    expected_r2 =  [-5200 NaN NaN NaN;-3500 NaN NaN NaN; 4637 NaN NaN NaN]
    expected_v2 = [5.200000e+00 NaN NaN NaN;6.400000e+00 NaN NaN NaN; 2.500000e+00 NaN NaN NaN]
    @test isapprox(r2, expected_r2, rtol=1e-6, atol=1e-8) # Compare calculated and expected position vectors
    @test isapprox(v2, expected_v2, rtol=1e-6, atol=1e-8) # Compare calculated and expected velocity vectors
    # 3. J2 perturbation case
    flag_drag = false
    flag_J2 = true
    r3, v3 = OrbitalPerturbations.cowell(r0, v0, Bc, t, flag_drag, flag_J2)
    # Define expected results for the test case
    expected_r3 = [-5200 5.871090e+03 6.068004e+03 6.192462e+03;
                    -3500 2.667597e+03 -7.452999e+02 -2.070397e+03;
                    4637 -1.232726e+04 -2.093067e+04 -2.292971e+04]
    expected_v3 = [5.200000e+00 -9.756580e-02 -4.993729e-01 -4.996058e-01;
                    6.400000e+00 -2.612848e+00 -2.423831e+00 -2.268180e+00; 
                    2.500000e+00 -4.739604e+00 -1.561818e+00 2.188937e-01]
    @test isapprox(r3, expected_r3, rtol=1e-6, atol=1e-8) # Compare calculated and expected position vectors
    @test isapprox(v3, expected_v3, rtol=1e-6, atol=1e-8) # Compare calculated and expected velocity vectors
    # 4. Drag and J2 perturbation case
    flag_drag = true
    flag_J2 = true
    r4, v4 = OrbitalPerturbations.cowell(r0, v0, Bc, t, flag_drag, flag_J2)
    # Define expected results for the test case
    expected_r4 = [-5200 NaN NaN NaN;-3500 NaN NaN NaN; 4637 NaN NaN NaN]
    expected_v4 = [5.200000e+00 NaN NaN NaN;6.400000e+00 NaN NaN NaN; 2.500000e+00 NaN NaN NaN]
    @test isapprox(r4, expected_r4, rtol=1e-6, atol=1e-8) # Compare calculated and expected position vectors
    @test isapprox(v4, expected_v4, rtol=1e-6, atol=1e-8) # Compare calculated and expected velocity vectors
end

# Run the tests
@testset "OrbitalPerturbations tests" begin 
    test_J2_acceleration()
    test_drag_acceleration()
    test_cowell()
end