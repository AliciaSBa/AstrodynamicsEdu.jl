# OrbitalPerturbations file

include("linearAlgebraTypes.jl")
include("idealTwoBodyProblem.jl")
include("astroConstants.jl")
using LinearAlgebra
using DifferentialEquations
using OrdinaryDiffEq

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
function cowell(r0, v0, bc, t, flag_drag, flag_J2)
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
    function propagate!(u, p)
        # Extract the parameters
        bc = p[1]
        flag_drag = p[2]
        flag_J2 = p[3]
    
        rr = u[1:3]
        vv = u[4:6]
        r = norm(rr)
        zeta = norm(rr) - R_Earth
    
        # Check drag
        if flag_drag==true
            acc_drag = drag_acceleration(zeta, vv, bc)
        else
            acc_drag = 0.0
        end
    
        # Check J2
        if flag_J2==true
            acc_J2 = J2_acceleration(rr)
        else
            acc_J2 = 0.0
        end
    
        # Compute RHS
        du = zeros(6)
        du[1:3] = vv
        du[4:6] = -mu_Earth / r^3 * rr .+ acc_drag .+ acc_J2
        #println("du: ", du)
        return du
    end

    p = [bc, flag_drag, flag_J2]

    # Define the problem
    odefun(u,p,t) = propagate!(u, p)

    # Integration
    u0 = [r0; v0]
    tspan = (t[1], t[end])
    prob = ODEProblem(odefun, u0, tspan, p)
    sol = solve(prob, Tsit5(), reltol = 1e-9, abstol = 1e-9)

    #r = permutedims(sol(t)[:, 1:3])
    #v = permutedims(sol(t)[:, 4:6])

    r = sol(t)[1:3,:]
    v = sol(t)[4:6,:]

    #println("r: ", r)
    #println("v: ", v)

    return r, v
end
