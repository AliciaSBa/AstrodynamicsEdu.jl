module OrbitalManeuvers

include("LinearAlgebraTypes.jl")
include("IdealTwoBodyProblem.jl")
using LinearAlgebra

export HohmannTransfer_circular, HohmannTransfer_elliptic, planeChange_circular, planeChange_apoapsis2apoapsis, 
        LowThrust_orbitRaising_circular, LowThrust_planeChange_circular, HohmannTransfer_tof

# Calculate the delta-v required to change the orbit of a spacecraft between two circular and coplanar orbits (Hohmann transfer)
function HohmannTransfer_circular(r1::Float64, r2::Float64, mu::Float64)
    # Input: 
        # r1: radius of the initial circular orbit
        # r2: radius of the final circular orbit
        # mu: gravitational parameter of the central body
    # Output:
        # deltaV: total delta-v required to change the orbit from r1 to r2
    # Calculate the semi-major axis of the transfer orbit
    at = (r1 + r2)/2
    # Calculate the delta-v required
    if r1 < r2
        # Calculate the delta-v1 at the first firing (at pericenter)
        v1 = circularVelocity(r1, mu)
        vtp = speed_visViva(r1, at, mu)
        deltaV1 = vtp - v1
        # Calculate the delta-v2 at the second firing (at apocenter)
        v2 = circularVelocity(r2, mu)
        vta = speed_visViva(r2, at, mu)
        deltaV2 = v2 - vta
    elseif r1 > r2
        # Calculate the delta-v1 at the first firing (at apocenter)
        v1 = circularVelocity(r1, mu)
        vta = speed_visViva(r1, at, mu)
        deltaV1 = v1 - vta
        # Calculate the delta-v2 at the second firing (at pericenter)
        v2 = circularVelocity(r2, mu)
        vtp = speed_visViva(r2, at, mu)
        deltaV2 = vtp - v2
    else
        deltaV1 = 0.0
        deltaV2 = 0.0
    end
    deltaV = deltaV1 + deltaV2
    return deltaV
end

# Calculate the delta-v required to change the orbit of a spacecraft between two elliptic and coplanar orbits (Hohmann transfer)
function HohmannTransfer_elliptic(rp1::Float64, e1::Float64, rp2::Float64, e2::Float64, mu::Float64)
    # Input: 
        # rp1: radius of the initial elliptic orbit at pericenter
        # e1: eccentricity of the initial elliptic orbit
        # rp2: radius of the final elliptic orbit at pericenter
        # e2: eccentricity of the final elliptic orbit
        # mu: gravitational parameter of the central body
    # Output:
        # deltaV: total delta-v required to change the orbit from rp1 to rp2
    # Calculate the semi-major axis of the initial orbit
    a1 = rp1/(1 - e1)
    # Calculate the semi-major axis of the final orbit
    a2 = rp2/(1 - e2)
        # Calculate the delta-v required
    if rp1 < rp2
        # Calculate the semi-major axis of the transfer orbit
        ra2 = a2*(1 + e2)
        at = (ra2 + rp1)/2
        # Calculate the delta-v1 at the first firing (at pericenter)
        v1 = speed_visViva(rp1, a1, mu)
        vtp = speed_visViva(rp1, at, mu)
        deltaV1 = vtp - v1
        # Calculate the delta-v2 at the second firing (at apocenter)
        vta = speed_visViva(ra2, at, mu)
        v2 = speed_visViva(ra2, a2, mu)
        deltaV2 = v2 - vta
    elseif rp1 >= rp2
        # Calculate the semi-major axis of the transfer orbit
        ra1 = rp1*(1-e1)/(1 + e1)
        at = (ra1 + rp2)/2
        # Calculate the delta-v1 at the first firing (at apocenter)
        v1 = speed_visViva(ra1, a1, mu)
        vta = speed_visViva(ra1, at, mu)
        deltaV1 = v1 - vta
        # Calculate the delta-v2 at the second firing (at pericenter)
        vtp = speed_visViva(rp2, at, mu)
        v2 = speed_visViva(rp2, a2, mu)
        deltaV2 = vtp - v2
    #else
        #deltaV1 = 0.0
        #deltaV2 = 0.0
    end
    deltaV = deltaV1 + deltaV2
    return deltaV
end

# Calculate the minimum delta-v required to change the inclination between two orbits (Plane Change Maneuver)
function planeChange_apoapsis2apoapsis(stateVector1::MyStateVector, stateVector2::MyStateVector, delta_i::Float64, mu::Float64)
    # Input: 
        # stateVector1: state vector of the initial orbit
        # stateVector2: state vector of the final orbit
        # delta_i: inclination change [rad]
        # mu: gravitational parameter of the central body
    # Output:
        # deltaV: total delta-v required to change the inclination from i1 to i2
    # Calculate the delta-v required
    v1 = velocityApoapsis(stateVector1, mu)
    v2 = velocityApoapsis(stateVector2, mu)
    deltaV = 2*sqrt(v1^2 + v2^2 - 2*v1*v2*cos(delta_i))
    return deltaV
end

# Calculate the minimum delta-v required to change the inclination between two orbits (Plane Change Maneuver)
function planeChange_apoapsis2apoapsis(coe1::MyCOE,coe2::MyCOE,mu::Float64)
    # Input: 
        # coe1: classical orbital elements of the initial orbit
        # coe2: classical orbital elements of the final orbit
        # mu: gravitational parameter of the central body
    # Output:
        # deltaV: total delta-v required to change the inclination from i1 to i2
    # Calculate the velocity at the apoapsis
    v1 = sqrt(mu/coe1.a*(1-coe1.e)/(1+coe1.e))
    v2 = sqrt(mu/coe2.a*(1-coe2.e)/(1+coe2.e))
    # Calculate the delta-v required
    deltaV = 2*sqrt(v1^2 + v2^2 - 2*v1*v2*cos(coe1.inc - coe2.inc))
    return deltaV
end

# Calculate the delta-v required to change the inclination between two circular orbits with equal radius (Plane Change Maneuver)
function planeChange_circular(ra::Float64, delta_i::Float64, mu::Float64)
    # Input: 
        # ra: radius of the circular orbits
        # delta_i: inclination change [rad]
        # mu: gravitational parameter of the central body
    # Output:
        # deltaV: total delta-v required to change the inclination from i1 to i2
    # Calculate the delta-v required
    v = circularVelocity(ra, mu)
    deltaV = 2*v*sin((delta_i)/2)
    return deltaV
end

# Calculate the delta-v required to go from a circular orbit to another circular and coplanar orbit (Low Thrust Orbit raising)
function LowThrust_orbitRaising_circular(r1::Float64,r2::Float64,mu::Float64)
    # Input: 
        # r1: radius of the initial circular orbit
        # r2: radius of the final circular orbit
        # mu: gravitational parameter of the central body
    # Output:
        # deltaV: total delta-v required to change the orbit from r1 to r2
    # Calculate the delta-v required
    v1 = circularVelocity(r1, mu)
    v2 = circularVelocity(r2, mu)
    deltaV = abs(v2 - v1)
    return deltaV
end

# Calculate the delta-v required to change inclination between two circular orbits with equal radius (Low Thrust Plane Change Maneuver)
function LowThrust_planeChange_circular(ra::Float64,delta_i::Float64,mu::Float64)
    # Input: 
        # ra: radius of the circular orbits
        # delta_i: inclination change [rad]
        # mu: gravitational parameter of the central body
    # Output:
        # deltaV: total delta-v required to change the orbit from r1 to r2
    # Calculate the delta-v required
    v = circularVelocity(ra, mu)
    deltaV = sqrt(2*v^2*(1-cos(pi/2*delta_i)))
    return deltaV
end

# Given a desired time of flight and the Classical Orbital Elements of two orbits, calculate the delta-v vector required 
# to perform a Hohmann transfer, as well as the the velocity vectors at the beginning and at the end of the transfer
function HohmannTransfer_tof(coe1::MyCOE, coe2::MyCOE, tof::Float64, mu::Float64)
    # Input: 
        # coe1: classical orbital elements of the initial orbit
        # coe2: classical orbital elements of the final orbit
        # tof: desired time of flight
        # mu: gravitational parameter of the central body
    # Output:
        # deltaV: delta-v vector required to perform the transfer
        # vT1: velocity vector at the beginning of the transfer
        # vT2: velocity vector at the end of the transfer
    # Calculate the initial and final position vectors
    state1 = COE_to_stateVector(coe1, mu)
    state2 = COE_to_stateVector(coe2, mu)
    r1 = state1.r
    r2 = state2.r
    coeT1,coeT2 = Lambert_solve(r1,r2,tof,mu,k0)
    if coeT1.e > 1.0
        println("The Hohmann transfer orbit is not possible. A hyperbolic orbit is required.")
        return NaN, NaN, NaN
    elseif isapprox(coeT1.e,1.0,atol=1e-4)
        println("The Hohmann transfer orbit is not possible. A parabolic orbit is required.")
        return NaN, NaN, NaN
    else
        # Calculate the initial and final velocity vectors
        stateT1 = COE_to_stateVector(coeT1, mu)
        stateT2 = COE_to_stateVector(coeT2, mu)
        vT1 = stateT1.v
        vT2 = stateT2.v
        # Calculate the delta-v vector
        deltaV1 = vT1 - state1.v
        deltaV2 = vT2 - state2.v
        deltaV = deltaV1 + deltaV2
        return deltaV, vT1, vT2
    end
end


end # module OrbitalManeuvering