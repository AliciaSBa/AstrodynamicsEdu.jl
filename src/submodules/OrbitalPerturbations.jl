module OrbitalPerturbations

#include("LinearAlgebraTypes.jl")
#include("IdealTwoBodyProblem.jl")
#include("AstroConstants.jl")
using AstrodynamicsEdu.LinearAlgebraTypes
using AstrodynamicsEdu.IdealTwoBodyProblem
using AstrodynamicsEdu.AstroConstants
using LinearAlgebra
using DifferentialEquations
#using OrdinaryDiffEq

export J2_acceleration, drag_acceleration, cowell

# Calculate the perturbation acceleration due to the J2 term (Earth)
function J2_acceleration(r::Vector)
    # Input:
        # r: position vector in the inertial frame
    # Output:
        # ap_J2: perturbation acceleration due to the J2 term
    # For the Earth
    J2 = J2_Earth
    # Calculate the perturbation acceleration
    r_norm = norm(r)
    z = r[3]
    ur = r/r_norm
    uz = [0.0,0.0,1.0]
    a_p_J2 = 3/2*J2*mu_Earth*R_Earth^2/r_norm^4*((5*(z/r_norm)^2 -1)*ur - 2*z/r_norm*k0.c0)
    #a_p_J2 = 3/2*J2*mu_Earth*R_Earth^2/r_norm^4*((5*dot(ur,uz)^2 -1)*ur - 2*dot(ur,uz)*uz)
    return a_p_J2
end

# Calculate the perturbation acceleration due to the Drag (Earth)
function drag_acceleration(alt::Float64,v::Vector,Bc::Float64)
    # Input:
        # alt: altitude of the satellite
        # v: velocity vector in the inertial frame
        # Bc: ballistic coefficient
    # Output:
        # a_p_drag: perturbation acceleration due to the drag
    # Assuming a simple isothermal module
    H = 50.0 # km (scale height)
    rho_0 = 1.0e-8*10^9# kg/km^3 (density at reference height)
    alt_0 = 100.0 # km (reference height)
    # Calculate the density at the given altitude
    rho = rho_0*exp(-(alt-alt_0)/H)
    # Calculate the perturbation acceleration
    v_norm = norm(v)
    if alt < 100.0
        #a_p_drag = [NaN,NaN,NaN]
        a_p_drag = NaN
    else
        a_p_drag = -1/(2*Bc)*rho*v_norm*v
    end
    return a_p_drag
end

# Calculate the perturbation acceleration due to the Solar Radiation Pressure

# Calculate the perturbation acceleration due to the Third Body B


# Establish an orbits propagator using the Cowell method
function cowell(r0::Vector{Float64}, v0::Vector{Float64}, Bc::Float64, tvec::Vector, flag_drag::Bool, flag_J2::Bool)
    # Input:
    # r0: initial position vector in the inertial frame
    # v0: initial velocity vector in the inertial frame
    # Bc: ballistic coefficient
    # t: time vector
    # flag_drag: Boolean that activates the drag perturbation when true
    # flag_J2: Boolean that activates the J2 perturbation when true
    # Output:
    # r: position vector in the inertial frame
    # v: velocity vector in the inertial frame

    # Define the differential equation for orbit propagation
    function propagate!(u,Bc,flag_drag,flag_J2)
        r = u[1:3]
        v = u[4:6]
        r_norm = norm(r)

        # Acceleration due to gravity (Keplerian motion)
        a_gravity = -mu_Earth * r / r_norm^3
    
        # Perturbations
        if flag_drag
            # Acceleration due to drag
            ap_drag = drag_acceleration(r_norm-R_Earth, v, Bc)
        else
            ap_drag = [0.0, 0.0, 0.0]
        end
        
        if flag_J2
            # Acceleration due to J2 perturbation
            ap_J2 = J2_acceleration(r)
        else
            ap_J2 = [0.0, 0.0, 0.0]
        end

        # Total acceleration
        du = zeros(6)
        du[1:3] = v
        du[4:6] = a_gravity + ap_drag + ap_J2
        return du
    end

    # Define the initial state
    u0 = [r0; v0]

    # Define the parameters
    parameters = [Bc, flag_drag, flag_J2]

    #ode_fun!(du, u, p, t) = propagate!(u, Bc, flag_drag, flag_J2)
    ode_fun!(u,p,t) = propagate!(u, Bc, flag_drag, flag_J2)

    # Define the problem
    #prob = ODEProblem(ode_fun!, u0, (0.0, t[end]))

    tspan = first(tvec), last(tvec)
    prob = ODEProblem(ode_fun!, u0, tspan)

    # Solve the problem
    #sol = solve(prob, Tsit5(), saveat = tvec, abstol = 1e-9, reltol = 1e-9)
    #sol = solve(prob,DP5(),abstol=1e-9,reltol=1e-9)
    #sol = solve(prob,DP5(),saveat=tvec,abstol=1e-9,reltol=1e-9)
    sol = solve(prob, Tsit5(), saveat = tvec)
    
    # Extract the solution
    r = sol[1:3, :]
    v = sol[4:6, :]

    return r, v

end



end # module OrbitalPerturbations