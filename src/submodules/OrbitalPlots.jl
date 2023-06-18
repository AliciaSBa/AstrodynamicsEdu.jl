module OrbitalPlots

#include("LinearAlgebraTypes.jl")
#include("IdealTwoBodyProblem.jl")
#include("AstroConstants.jl")
using AstrodynamicsEdu.LinearAlgebraTypes
using AstrodynamicsEdu.IdealTwoBodyProblem
using AstrodynamicsEdu.AstroConstantss
using LinearAlgebra
using Plots
#using DifferentialEquations
#using Printf
#using Plots.Measures

export plot_orbit_perifocal, plot_orbit_ECI, plot_porkchop

################## WORKING PLOTS ######################

# Plot an elliptical orbit in the perifocal frame (2D)
function plot_orbit_perifocal_E(coe::MyCOE, planet_radius::Float64)
    theme(:juno)
    gr()
    # Extract the orbital elements
    e = coe.e
    p = coe.p
    theta = coe.theta

    # Calculate the position of the planet at the focus
    planet_x = 0.0
    planet_y = 0.0

    # Generate points on the orbit
    th = range(0, stop=2π, length=1000)
    r = p ./ (1 .+ e.*cos.(th))

    # Calculate x and y coordinates in the perifocal frame
    x = r.*cos.(th)
    y = r.*sin.(th)

    # Create the plot
    plot(x, y, aspect_ratio=:equal, linewidth=1.5, legend=false, linecolor=:aliceblue, grid=false, framestyle=:box,dpi=600)
    title!("Orbit in the Perifocal plane")
    xlabel!("X[km]")
    ylabel!("Y[km]")

    # Plot the planet at the focus
    #scatter!([planet_x], [planet_y], marker=:circle, color=:mediumaquamarine, markersize=(planet_radius/1000),markerstrokewidth = 0.0)
    xE(t) = planet_radius*cos(t)
    yE(t) = planet_radius*sin(t)
    plot!(xE, yE, 0, 2π, color=:mediumaquamarine, linewidth=2.0, fill=true, fillcolor=:mediumaquamarine)

    # Calculate the position of the satellite in the orbit
    r_sat = p / (1 + e*cos(theta))
    satellite_x = r_sat*cos(theta)
    satellite_y = r_sat*sin(theta)

    # Plot the satellite in the orbit
    scatter!([satellite_x], [satellite_y], marker=:circle, color=:sandybrown, markersize = 5.0, markerstrokewidth = 0.0)
end

# Plot a hyperbolic orbit in the perifocal frame (2D)
function plot_orbit_perifocal_H2(coe::MyCOE, planet_radius::Float64)
    theme(:juno)
    gr()
    # Extract the orbital elements
    a = coe.a
    e = coe.e

    # Calculate the semi-minor axis of the hyperbola
    b = a * sqrt(e^2 - 1)

    # Calculate the distance from the center to the focus
    c = sqrt(a^2 + b^2)

    # Determine which focus to use
    if e > 1
        focus_x = -c
    else
        focus_x = c
    end

    focus_y = 0.0

    # Calculate the asymptote angle
    theta_asymptote = acos(-1 / e)

    # Generate points on the orbit
    t = range(-theta_asymptote, stop=theta_asymptote, length=2000)
    r = a * (e^2 - 1) ./ (1 .+ e * cosh.(t))

    # Calculate x and y coordinates in the perifocal frame for the first branch of the hyperbola
    x1 = r.*cosh.(t)
    y1 = r.*sinh.(t)

    # Calculate x and y coordinates in the perifocal frame for the second branch of the hyperbola
    x2 = -x1
    y2 = -y1

    # Create the plot
    plot(x1, y1, aspect_ratio=:equal, linewidth=1.5, legend=false, linecolor=:aliceblue, grid=false, framestyle=:box, dpi=600)
    #plot!(x2, y2, linewidth=1.5, linecolor=:aliceblue)

    title!("Orbit in the Perifocal plane")
    xlabel!("X[km]")
    ylabel!("Y[km]")

    # Plot the planet at the focus
    # scatter!([focus_x], [focus_y], marker=:circle, color=:mediumaquamarine, markersize=(planet_radius/1000),markerstrokewidth = 0.0)
    xE(t) = planet_radius*cos(t)
    yE(t) = planet_radius*sin(t)
    plot!(xE, yE, 0, 2π, color=:mediumaquamarine, linewidth=2.0, fill=true, fillcolor=:mediumaquamarine)

    # Calculate the position of the satellite in the orbit
    r_theta = a * (e^2 - 1) ./ (1 .+ e * cosh.(theta))
    satellite_x = r_theta.*cosh.(theta)
    satellite_y = r_theta.*sinh.(theta)
    #satellite_x = a * (e - cosh(coe.theta))
    #satellite_y = a * sqrt(e^2 - 1) * sinh(coe.theta)

    # Plot the satellite in the orbit
    scatter!([satellite_x], [satellite_y], marker=:circle, color=:sandybrown, markersize = 5.0, markerstrokewidth = 0.0)
end

function plot_orbit_perifocal_H(coe::MyCOE, planet_radius::Float64)
    theme(:juno)
    gr()
    # Extract the orbital elements
    a = coe.a
    e = coe.e

    # Calculate the semi-minor axis of the hyperbola
    b = a * sqrt(e^2 - 1)

    # Calculate the distance from the center to the focus
    c = sqrt(a^2 + b^2)

    # Determine which focus to use
    planet_x = a*e-a
    planet_y = 0.0

    # Calculate the asymptote angle
    theta_asymptote = acos(-1 / e)

    # Generate points on the orbit
    t = range(-theta_asymptote, stop=theta_asymptote, length=2000)
    r = a * (e^2 - 1) ./ (1 .+ e * cosh.(t))

    # Calculate x and y coordinates in the perifocal frame
    x1 = abs(a) .* (e .- cosh.(t))
    y1 = abs(b) .* sinh.(t)
    #x2 = -x1 .- 2*a*e
    #y2 = y1

    # Create the plot
    plot(x1, y1, aspect_ratio=:equal, linewidth=1.5, legend=false, linecolor=:aliceblue, grid=false, framestyle=:box, dpi=600)
    #plot!(x2, y2, linewidth=1.5, linecolor=:aliceblue)

    title!("Orbit in the Perifocal plane")
    xlabel!("X[km]")
    ylabel!("Y[km]")

    # Plot the planet at the focus
    #scatter!([planet_x], [planet_y], marker=:circle, color=:mediumaquamarine, markersize=(planet_radius/1000), markerstrokewidth = 0.0)
    xE(t) = planet_radius*cos(t)
    yE(t) = planet_radius*sin(t)
    plot!(xE, yE, 0, 2π, color=:mediumaquamarine, linewidth=2.0, fill=true, fillcolor=:mediumaquamarine)

    # Calculate the position of the satellite in the orbit
    satellite_x = abs(a) .* (e .- cosh.(theta))
    satellite_y = abs(b) .* sinh.(theta)
    #satellite_x = a * (e - cosh(coe.theta))
    #satellite_y = a * sqrt(e^2 - 1) * sinh(coe.theta)

    # Plot the satellite in the orbit
    scatter!([satellite_x], [satellite_y], marker=:circle, color=:sandybrown, markersize = 5.0, markerstrokewidth = 0.0)
end

# Plot a parabolic orbit in the perifocal frame (2D)
function plot_orbit_perifocal_P(coe::MyCOE,planet_radius::Float64)
    theme(:juno)
    gr()
    # Extract the orbital elements
    p = coe.p
    e = coe.e
    theta = coe.theta

    # Calculate the trajectory of the Parabola
    t = range(-2*pi/3, stop=2*pi/3, length=2000)
    r = p ./ (1 .+ e * cos.(t))
    x = r.*cos.(t)
    y = r.*sin.(t)

    # Calculate the position of the satellite in the orbit
    r_theta = p ./ (1 .+ e * cos.(theta))
    satellite_x = r_theta.*cos.(theta)
    satellite_y = r_theta.*sin.(theta)

    # Calculate the position of the planet at the focus
    planet_x = 0.0
    planet_y = 0.0

    # Create the plot
    plot(x, y, aspect_ratio=:equal, linewidth=1.5, legend=false, linecolor=:aliceblue, grid=false, framestyle=:box, dpi=600)
    title!("Orbit in the Perifocal plane")
    xlabel!("X[km]")
    ylabel!("Y[km]")

    # Plot the planet at the focus
    scatter!([planet_x], [planet_y], marker=:circle, color=:mediumaquamarine, markersize=(planet_radius/1000), markerstrokewidth = 0.0)
    xE(t) = planet_radius*cos(t)
    yE(t) = planet_radius*sin(t)
    plot!(xE, yE, 0, 2π, color=:mediumaquamarine, linewidth=2.0, fill=true, fillcolor=:mediumaquamarine)

    # Plot the satellite in the orbit
    scatter!([satellite_x], [satellite_y], marker=:circle, color=:sandybrown, markersize = 5.0, markerstrokewidth = 0.0)
end

function plot_orbit_perifocal(coe::MyCOE,planet_radius::Float64)
    if isapprox(coe.e, 1.0, atol=1e-4)
        plot_orbit_perifocal_P(coe, planet_radius)
    elseif coe.e > 1.0
        plot_orbit_perifocal_H(coe, planet_radius)
    elseif coe.e < 1.0 && coe.e >= 0.0
        plot_orbit_perifocal_E(coe, planet_radius)
    else
        error("Invalid eccentricity")
    end
end


function plot_orbit_ECI(coe::MyCOE,planet_radius::Float64, mu::Float64)
    theme(:juno)
    plotly()
    #gr()

    function sphere(radius::Float64, resolution::Int)
        theta = range(0, stop = 2π, length = resolution)
        phi = range(0, stop = π, length = resolution)
        x = zeros(resolution, resolution)
        y = zeros(resolution, resolution)
        z = zeros(resolution, resolution)
        for i=1:resolution
            for j=1:resolution
                x[i, j] = radius * sin(phi[i]) * cos(theta[j])
                y[i, j] = radius * sin(phi[i]) * sin(theta[j])
                z[i, j] = radius * cos(phi[i])
            end
        end

        surface(x, y, z, color = :mediumaquamarine, colorbar=:none, legend=:none, framestyle=:box, dpi=600)
    end
    
    # Plot the planet
    sphere(planet_radius, 50)

    # Plot the elliptical orbit
    r3D = zeros(3, 1000)
    theta_it = range(0, stop = 2π, length = 1000)
    for i = 1:1000
        theta = theta_it[i]
        coe_it = MyCOE(RAAN=coe.RAAN, inc=coe.inc, omega=coe.omega, theta=theta, p=coe.p, e=coe.e)
        state_it = COE_to_stateVector(coe_it,mu)
        r3D[:, i] = state_it.r.c0[1:3]
    end

    plot3d!(r3D[1, :], r3D[2, :], r3D[3, :], aspect_ratio = 1, linecolor=:aliceblue,linewidth = 1.0)

    # Plot the satellite
    state_Sat = COE_to_stateVector(coe,mu)
    r_Sat = state_Sat.r.c0
    scatter3d!([r_Sat[1]], [r_Sat[2]], [r_Sat[3]], marker=:circle, color=:sandybrown, markersize = 2.5, markerstrokealpha=0.0)
    title!("Orbit in the ECI frame")
    xlabel!("X[km]")
    ylabel!("Y[km]")
    zlabel!("Z[km]")


end

################## IN DEVELOPMENT ######################

function plot_orbit_ECI_H(coe::MyCOE,planet_radius::Float64, mu::Float64)
    theme(:juno)
    plotly()
    #gr()

    function sphere(radius::Float64, resolution::Int)
        theta = range(0, stop = 2π, length = resolution)
        phi = range(0, stop = π, length = resolution)
        x = zeros(resolution, resolution)
        y = zeros(resolution, resolution)
        z = zeros(resolution, resolution)
        for i=1:resolution
            for j=1:resolution
                x[i, j] = radius * sin(phi[i]) * cos(theta[j])
                y[i, j] = radius * sin(phi[i]) * sin(theta[j])
                z[i, j] = radius * cos(phi[i])
            end
        end

        surface(x, y, z, color = :mediumaquamarine, colorbar=:none, legend=:none, framestyle=:box, dpi=600)
    end
    
    # Plot the planet
    sphere(planet_radius, 50)

    # Plot the hyperbolic orbit
    r3D = zeros(3, 1000)
    theta_asymptote = acos(-1/coe.e)
    theta_it = range(-theta_asymptote, stop=theta_asymptote, length=1000)
    #theta_it = range(0, stop = 2π, length = 1000)
    for i = 1:1000
        theta = theta_it[i]
        coe_it = MyCOE(RAAN=coe.RAAN, inc=coe.inc, omega=coe.omega, theta=theta, p=coe.p, e=coe.e)
        state_it = COE_to_stateVector(coe_it,mu)
        r3D[:, i] = state_it.r.c0[1:3]
    end

    plot3d!(r3D[1, :], r3D[2, :], r3D[3, :], aspect_ratio = 1, linecolor=:aliceblue,linewidth = 1.0)


    # Plot the satellite
    state_Sat = COE_to_stateVector(coe,mu)
    r_Sat = state_Sat.r.c0
    scatter3d!([r_Sat[1]], [r_Sat[2]], [r_Sat[3]], marker=:circle, color=:sandybrown, markersize = 2.5, markerstrokealpha=0.0)
    title!("Orbit in the ECI frame")
    xlabel!("X[km]")
    ylabel!("Y[km]")
    zlabel!("Z[km]")


end

function plot_orbit_ECI_P(coe::MyCOE,planet_radius::Float64, mu::Float64)
    theme(:juno)
    plotly()
    #gr()

    function sphere(radius::Float64, resolution::Int)
        theta = range(0, stop = 2π, length = resolution)
        phi = range(0, stop = π, length = resolution)
        x = zeros(resolution, resolution)
        y = zeros(resolution, resolution)
        z = zeros(resolution, resolution)
        for i=1:resolution
            for j=1:resolution
                x[i, j] = radius * sin(phi[i]) * cos(theta[j])
                y[i, j] = radius * sin(phi[i]) * sin(theta[j])
                z[i, j] = radius * cos(phi[i])
            end
        end

        surface(x, y, z, color = :mediumaquamarine, colorbar=:none, legend=:none, framestyle=:box, dpi=600)
    end
    
    # Plot the planet
    sphere(planet_radius, 50)

    # Plot the parabolic orbit
    r3D = zeros(3, 1000)
    theta_it = range(-2*pi/3, stop=2*pi/3, length = 1000)
    for i = 1:1000
        theta = theta_it[i]
        coe_it = MyCOE(RAAN=coe.RAAN, inc=coe.inc, omega=coe.omega, theta=theta, p=coe.p, e=coe.e)
        state_it = COE_to_stateVector(coe_it,mu)
        r3D[:, i] = state_it.r.c0[1:3]
    end

    plot3d!(r3D[1, :], r3D[2, :], r3D[3, :], aspect_ratio = 1, linecolor=:aliceblue,linewidth = 1.0)

    # Plot the satellite
    state_Sat = COE_to_stateVector(coe,mu)
    r_Sat = state_Sat.r.c0
    scatter3d!([r_Sat[1]], [r_Sat[2]], [r_Sat[3]], marker=:circle, color=:sandybrown, markersize = 2.5, markerstrokealpha=0.0)
    title!("Orbit in the ECI frame")
    xlabel!("X[km]")
    ylabel!("Y[km]")
    zlabel!("Z[km]")

end

function plot_porkchop(t_launch,t_travel,pork_plot, t_synodic)
    # Input:
        # t_launch: vector of launch dates [days]
        # t_travel: vector of travel times [days]
        # pork_plot: matrix of C3 values [km/s] (length(t_launch) x length(t_travel))
        # t_synodic: synodic period [days]

    # Find the minimum energy and get the launch time and travel time
    minC3 = minimum(skipmissing(isnan(x) ? missing : x for x in pork_plot))
    positions = findall(x -> x == minC3, pork_plot)

    t_launch_min = t_launch[positions[1][2]]
    t_travel_min = t_travel[ positions[1][1]] # days

    # Calculate the synodic period to determine the next optimal launch window
    t_launch_next = t_launch_min + t_synodic
    t_travel_next = t_travel_min
    
    theme(:lime)
    gr(size =(2000,1000))

    # Plot the porkchop
    figurePlot = contourf(t_launch,t_travel,pork_plot, clims=(minC3,10*minC3))
    scatter!([t_launch_min],[t_travel_min],marker=true,markersize=5,markercolor="red",markerstrokecolor="red",markerstrokewidth=0.0,label="Minimum Energy")
    scatter!([t_launch_next],[t_travel_next],marker=true,markersize=5,markercolor="red",markerstrokecolor="red",markerstrokewidth=0.0,label="Next Optimal Launch Window")
    title!("C3 Porkchop Plot [km^2/s^2]")
    xlabel!("Launch time [days]")
    ylabel!("Travel time [days]")
    zlabel!("C3 [km/s]")
end



end # module OrbitalPlots


