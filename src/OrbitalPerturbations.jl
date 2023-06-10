module OrbitalPerturbations

using LinearAlgebra
using LinearAlgebraTypes
using IdealTwoBodyProblem
using DifferentialEquations
using AstroConstants

export J2_acceleration, drag_acceleration, cowell_propagator


# Calculate the perturbation acceleration due to the J2 term (Earth)
function J2_acceleration(r::MyVector)
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
    a_p_J2 = 3/2*J2*mu_Earth*R_Earth^2/r_norm^4*((5*(z/r_norm)^2 -1)*ur - 2*z/r_norm*k0)
    return a_p_J2
end

# Calculate the perturbation acceleration due to the Drag (Earth)
function drag_acceleration(alt::Float64,v::MyVector,Bc::Float64)
    # Input:
        # alt: altitude of the satellite
        # v: velocity vector in the inertial frame
        # Bc: ballistic coefficient
    # Output:
        # a_p_drag: perturbation acceleration due to the drag
    # Assuming a simple isothermal module
    H = 50.0 # km (scale height)
    rho_0 = 1.0e-8 # kg/m^3 (density at reference height)
    alt_0 = 100.0 # km (reference height)
    # Calculate the density at the given altitude
    rho = rho_0*exp(-(alt-alt_0)/H)
    # Calculate the perturbation acceleration
    v_norm = norm(v)
    if alt < 100.0
        a_p_drag = MyVector(NaN,NaN,NaN)
    else
        a_p_drag = 1/(2*Bc)*rho*v_norm*v
    end
    return a_p_drag
end

# Calculate the perturbation acceleration due to the Solar Radiation Pressure

# Calculate the perturbation acceleration due to the Third Body B


# Establish an orbits propagator using the Cowell method
function cowell_propagator(r0::Vector{Float64}, v0::Vector{Float64}, Bc::Float64, t::Vector, flag_drag::Bool, flag_J2::Bool)
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
    function orbit_eqn!(du, u, p, t)
        r = u[1:3]
        v = u[4:6]
        r_norm = norm(r)
        
        # Acceleration due to gravity (Keplerian motion)
        a_gravity = -mu_Earth * u[1:3] / r_norm^3
        
        # Perturbations
        if flag_drag
            # Acceleration due to drag
            ap_drag = drag_acceleration(norm(r0) - R_Earth, MyVector(v0), Bc)
        else
            ap_drag = MyVector(0.0, 0.0, 0.0)
        end
        
        if flag_J2
            # J2 perturbation
            ap_J2 = J2_acceleration(MyVector(r0))
        else
            ap_J2 = MyVector(0.0, 0.0, 0.0)
        end
        
        # Total acceleration
        du[1:3] = v
        du[4:6] = a_gravity + ap_drag.c0 + ap_J2.c0
    end
    
    # Define the initial state
    u0 = [r0; v0]

    # Define the problem
    prob = ODEProblem(orbit_eqn!, u0, (0.0, t[end]))

    # Solve the problem
    sol = solve(prob, Tsit5(), saveat = t)
    
    # Extract the solution
    r = sol[1:3, :]
    v = sol[4:6, :]
    
    return r, v
end




end # module OrbitalPerturbations