# Orbital Visualization

One of the main objectives of this library is not only to allow the computation of orbital dynamics but also the visualization of orbits in an appealing and engaging way. 

The orbitalPlots.jl file provides a powerful toolset for generating high-quality orbit representations using the Plots package in Julia. By leveraging the capabilities of Plots, this module enables the creation of visually appealing and informative plots of various orbital elements and parameters

Upon initial use, it is worth noting that the orbitalPlots.jl file, and therefore the package, may take some time to load due to its dependence on the heavy Plots package. However, once loaded, the module functions smoothly and efficiently. In some cases, when writing code outside of the terminal, it may be necessary to include the `display()` function to ensure that the plots are rendered and displayed correctly.

Because of the relevance and regularity of these reference frames, this module allows to plot orbits in both the ECI and the perifocal reference frames. In the perifocal reference frame, the plots are done in the perifocal plane, in 2D. It uses the `gr()` backend to visualize them. Whereas, the plots in ECI reference frame are in 3D, and use the `plotly()`. The reason behind using this backend for the ECI plots instead of `gr()`, is that `plotly()` allows the user not only to visualize but also to interact with it. One is able to move around the plot and observe the orbit from a different perspective. For aesthetic reasons, both use `:juno` theme from PlotThemes.

The following examples showcase the visualization of orbits in the perifocal (PF) coordinate system, providing insights into the shape, size, and orientation of the orbits. 

![PFplotE](https://user-images.githubusercontent.com/115453770/252145648-41c6bdc4-0f6c-4f74-a0f5-28fb103ee045.png)
![PFplotP](https://user-images.githubusercontent.com/115453770/252145662-0bd3a9d6-ca44-44a9-a334-cd0974eae2a8.png)
![PFplotH](https://user-images.githubusercontent.com/115453770/252145661-eb68951a-56eb-43f5-8e64-484efaee26cb.png)


Additionally, an example in the Earth-Centered Inertial (ECI) coordinate system is also presented in the following figure. It shows a satellite in a LEO (Low Earth Orbit) around Earth. This coordinate system is commonly used in space missions and satellite tracking.

![ECIplotE](https://user-images.githubusercontent.com/115453770/252146523-2eb3a46c-9501-4035-930a-5f13375b5aee.png)

These examples serve as valuable references for understanding the behavior and characteristics of various types of orbits. Whether it's a circular orbit, elliptical orbit, or more complex trajectories, the orbitalPlots.jl file empowers students to gain deeper insights into orbital dynamics through visually compelling representations.

In the following table it can be seen the functions needed to do each plot, as well as the inputs required. There also appears a function to obtain a characteristic energy (C3) porkchop plot called `plot_porkchop`, an example of its usage and how it looks is shown in the Example Problem 2.

| **Function Name**          | **Inputs**                                            | **Functionality**                                                                                               |
|------------------------|---------------------------------------------------|-------------------------------------------------------------------------------------------------------------|
| plot\_orbit\_perifocal   | coe::MyCOE  planet\_radius::Float64              | Plot the orbit around a planet and the position of the satellite on it in the perifocal plane (2D)        |
| plot\_orbit\_ECI         | coe::MyCOE  planet\_radius::Float64  mu::Float64 | Plot the orbit around a planet and the position of the satellite on it in the ECI frame (3D)              |
| plot\_porkchop          | t\_launch::Vector [days]  t_travel::Vector [days]  pork\_plot::Matrix [km^2/s^2]  t\_synodic::Number [days] | Plot the C3 porkchop from a matrix of length(t\_launch) x length(t\_travel) containing the characteristic energy |
