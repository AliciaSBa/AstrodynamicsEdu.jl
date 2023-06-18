# ORBITAL MANEUVERS TESTS

#=
using LinearAlgebra
include("linearAlgebraTypes.jl")
include("idealTwoBodyProblem.jl")
include("orbitalManeuvers.jl")
inlcude("astroConstants.jl")
using Test
=#

# Test the functions in the Orbital Maneuvers module

# Homework 8. Case II. Hohmann transfer between two circular orbits from LEO to GEO at an inclination of 28.5 deg
function test_HohmannCircular()
    #Data:
    rA = R_Earth + 200.0 # km (LEO)
    rB = R_Earth + 35786.0 # km (GEO)
    deltaV = HohmannTransfer_circular(rA, rB, mu_Earth)
    @test isapprox(deltaV,3.94,rtol=1e-2, atol=1e-4) # km/s
end

# Homework 8. Case II. A change of plane is performed to change the inclination of a GEO orbit from 28.5 deg to 0 deg
function test_PlaneChangeCircular()
    # Data:
    rB = R_Earth + 35786.0 # km (GEO)
    delta_i = deg2rad(28.5)
    deltaV = planeChange_circular(rB, delta_i, mu_Earth)
    @test isapprox(deltaV,1.52,rtol=1e-2, atol=1e-4) # km/s
end

# Homework 7. Problem 1.a. Hohmann transfer between elliptic orbits. (adapted from Curtis)
function test_HohmannElliptic()
    # Data:
    rp1 = 7000.0 # km
    e1 = 0.2 
    rp2 = 37000.0 # km
    e2 = 0.5
    # Calculate the minimum delta-v required to change the orbit from rp1 to rp2
    deltaV = HohmannTransfer_elliptic(rp1, e1, rp2, e2, mu_Earth)
    @test isapprox(deltaV,2.77,rtol=1e-3, atol=1e-4) # km/s
end

# EXAMPLE 3. Plane change between two orbits
function test_planeChangeElliptic()
    # Data:
    # Case 1: Using COEs
    coe1 = MyCOE(RAAN=0.0, inc=0.0, omega=0.0, a=7000.0, e=0.2, theta=0.0)
    coe2 = MyCOE(RAAN=0.0, inc=deg2rad(25.7), omega=0.0, a=8000.0, e=0.4, theta=0.0)
    deltaVc = planeChange_apoapsis2apoapsis(coe1, coe2, mu_Earth)
    @test isapprox(deltaVc,5.66,rtol=1e-3, atol=1e-4) # km/s
    # Case 2: Using state vectors
    state1 = COE_to_stateVector(coe1, mu_Earth)
    state2 = COE_to_stateVector(coe2, mu_Earth)
    delta_i = coe2.inc - coe1.inc # rad
    deltaVs = planeChange_apoapsis2apoapsis(state1, state2, delta_i, mu_Earth)
    @test isapprox(deltaVs,deltaVc,rtol=1e-3, atol=1e-4) # km/s
end

# EXAMPLE 4: Low Thrust Orbit raising between 2 circular and coplanar orbits
function test_LowThrustOrbitRaising()
    #Data:
    rA = R_Earth + 200.0 # km (LEO)
    rB = R_Earth + 35786.0 # km (GEO)
    deltaV = LowThrust_orbitRaising_circular(rA, rB, mu_Earth)
    @test isapprox(deltaV,4.71,rtol=1e-2, atol=1e-4) # km/s
end

# EXAMPLE 5: Low Thrust Plane Change Maneuver between two circular orbits with equal radius
function test_LowThrustPlaneChange()
    # Data:
    ra = R_Earth + 10000.0 # km (GEO)
    delta_i = deg2rad(25.0)
    deltaV = LowThrust_planeChange_circular(ra, delta_i, mu_Earth)
    @test isapprox(deltaV,3.32,rtol=1e-2, atol=1e-4) # km/s
end

# EXAMPLE 6: Given a time of flight and the initial and final orbits, calculate the delta-v required to perform a Hohmann transfer
function test_HohmannTransfer_tof()
    coe1 = MyCOE(RAAN=0.0, inc=0.0, omega=0.0, a=7000.0, e=0.2, theta=0.0)
    coe2 = MyCOE(RAAN=0.0, inc=0.0, omega=0.0, a=8000.0, e=0.4, theta=3.14)
    mu = 398600.0 # km^3/s^2
    tof = 36000.0 # s
    deltaV,v1,v2 = HohmannTransfer_tof(coe1, coe2, tof, mu)
    @test isapprox(norm(deltaV),11.16,rtol=1e-3, atol=1e-4) # km/s
    @test isapprox(v1,MyVector([5.58,9.74,0.0]),rtol=1e-3, atol=1e-4) # km/s
    @test isapprox(v2,MyVector([5.57,-4.88,0.0]),rtol=1e-3, atol=1e-4) # km/s
end


# Run the tests
@testset "OrbitalManeuvers tests" begin 
    test_HohmannCircular()
    test_HohmannElliptic()
    test_PlaneChangeCircular()
    test_planeChangeElliptic()
    test_LowThrustOrbitRaising()
    test_LowThrustPlaneChange()
    test_HohmannTransfer_tof()
end





