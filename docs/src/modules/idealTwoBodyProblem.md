# Ideal Two Body Problem

The study of the ideal two-body problem forms the core of our codebase, serving as the foundation for more complex algorithms and calculations in orbital dynamics. This module introduces two key types: the state vector and the coe (classical orbital elements); and provides a comprehensive set of functions that enable precise analysis and prediction of celestial motion.

The ideal two-body problem focuses on the motion of two point masses, assuming no external influences from other celestial bodies. Although this simplification may appear trivial, it serves as a fundamental building block for understanding and modeling more intricate orbital dynamics scenarios.

Within the idealTwoBodyProblem.jl file, a range of functions has been developed to address various aspects of the ideal two-body problem. These functions form the backbone of orbital computations, providing essential tools for the calculation of critical parameters and transformations.

The figure below shows the two new data types defined in the ideal two-body problem file: MyStateVector and MyCOE. Below each appears their store objects and their type (indicated as `fieldName::Type`). These two data types are crucial for astrodynamics, as they allow to completely describe the motion of a body in space.

![I2BPTypes](https://user-images.githubusercontent.com/115453770/252145528-b9c0a030-89e4-4cf8-a899-0133fc8ad94c.png)

It is important to understand what each of these types entails and how to construct them:

1. **MyStateVector**: stores the position vector (`r`) and the velocity vector (`v`) of a body in space. Although there is no general imposition on the units, because of the nature of distances in space, the position of the satellite/spacecraft is usually given in km, while the velocity is given in km/s. Providing the coordinates in these units ensures the correct functioning of all the functions in the code. There is only one possible way to construct a MyStateVector object, and it is as follows:
   - `MyStateVector(r, v)`, where both `r` and `v` must be MyVector objects.

2. **MyCOE**: stores all six classical orbital elements (`Ω, i, ω, a, e, θ`), and `p`, so that the orbit can always be well defined. This is because of the peculiarity of the parabolic orbit, where `a = ∞`. However, for the construction, only six of the elements are given, where the user can choose to either provide the semi-major axis `a` or the semi-latus rectum `p`, but never both. However, this seventh element can always be accessed as it is computed internally. It is very important to note that all the angles must be given in radians, which is not optional. Additionally, all the parameters must be of type Float64. Examples for each constructor are provided below:
   - `MyCOE(RAAN = 0.18, inc = 0.26, omega = 0.30, a = 12000.0, e = 0.54, theta = π)`, where `a` is provided and `p` is computed internally.
   - `MyCOE(RAAN = 0.52, inc = 0.0, omega = 0.92, p = 8000.0, e = 1.0, theta = 2.09)`, where `p` is provided and `a` is computed internally.

This first table presents an overview of the functions available within the ideal two-body problem file. These functions encompass a wide range of functionalities, including the computation of the angular momentum vector, the conversion between state vectors and classical orbital elements, and the propagation of state vectors over time. With a total of 26 functions, this comprehensive set ensures that users have the necessary tools to accurately analyze and simulate the motion of celestial objects.

| **Function Name** | **Inputs** | **Outputs** | **Functionality** |
|-------------------|------------|-------------|-------------------|
| angularMomentumVector | stateVector::MyStateVector | h::MyVector | Obtain `h` vector |
| mechanicalEnergy | stateVector::MyStateVector, mu::Float64 | xi::Float64 | Obtain `ξ` |
| eccentricityVector | stateVector::MyStateVector, mu::Float64 | e::MyVector | Obtain `e` vector |
| semilatusRectum | stateVector::MyStateVector, mu::Float64 | p::Float64 | Obtain `p` |
| periodOrbit | stateVector::MyStateVector, mu::Float64, OrbitType::String | period::Float64 | Obtain period `τ` [s] |
| perifocalBasis | stateVector::MyStateVector, mu::Float64 | perifocalBasis::MyBasis | Perifocal Basis |
| periapsis | stateVector::MyStateVector, mu::Float64 | r\_p::Float64 | Obtain `r_p` |
| apoapsis | stateVector::MyStateVector, mu::Float64 | r\_a::Float64 | Obtain `r_a` |
| trajectoryEquation | stateVector::MyStateVector, mu::Float64, theta::Float64 | trajectoryEquation::Float64 | Obtain `r` |
| inclination | stateVector::MyStateVector | inc::Float64 | Obtain `i` |
| semiMajorAxis | stateVector::MyStateVector, mu::Float64 | a::Float64 | Obtain `a` |
| orbitType | stateVector::MyStateVector, mu::Float64 | orbitType::String | Obtain orbit type: - Elliptic - Hyperbolic - Parabolic - Circular |
| trueAnomaly | stateVector::MyStateVector, mu::Float64 | theta::Float64 | Obtain `θ` [rad] |
| stateVector\_to\_COE | stateVector::MyStateVector, mu::Float64 | coe::MyCOE | `(r, v) → COE` |
| COE\_to\_stateVector | coe::MyCOE, mu::Float64 | stateVector::MyStateVector | `COE → (r, v)` |
| escapeVelocity | r::MyVector, mu::Float64 | v\_esc::Float64 | Obtain `v_e` |
| circularVelocity | 2 OPTIONS: - r::MyVector, mu::Float64 - r\_norm::Float64, mu::Float64 | v\_circ::Float64 | Obtain `v_c` |
| velocityPeriapsis | stateVector::MyStateVector, mu::Float64 | v_p::Float64 | Obtain `v_p` |
| velocityApoapsis | stateVector::MyStateVector, mu::Float64 | v_a::Float64 | Obtain `v_a` |
| velocityHyperbolicAsymptote | stateVector::MyStateVector, mu::Float64 | v_h::Float64 | Obtain `v_h` |
| speed\_visViva | r\_norm::Float64, a::Float64, mu::Float64 | v::Float64 | Obtain `v` |
| rightAscensionOfTheAscendingNode | stateVector::MyStateVector | RAAN::Float64 | Obtain `Ω` |
| argumentOfPeriapsis | stateVector::MyStateVector, mu::Float64 | omega::Float64 | Obtain `ω` |
| vector\_in\_ECI\_to\_perifocal | ECI\_vector::Vector, coe::MyCOE | perifocalVector::MyVector | `r_ECI → r_PQW` |
| vector\_in\_perifocal\_to\_ECI | perifocalVector::Vector, coe::MyCOE | ECI\_vector::MyVector | `r_PQW → r_ECI` |
| propagate\_StateVector | stateVector::MyStateVector, mu::Float64, delta\_t::Float64 | stateVector2::MyStateVector | Propagate `(r, v)` by `Δt` |


Additionally, as displayed on the second table, the file incorporates functions specific to Kepler's problem, which focuses on the motion of bodies under the influence of a central gravitational force. The aim of this part is to be able to relate time and position. In order to achieve that, it is necessary to be able to convert between the different anomalies at ease. These conversions differ depending on the orbit type (elliptic, parabolic or hyperbolic), and obtaining them can become quite tricky. Because for the path from the time to the auxiliary anomaly, an iterative numerical method is needed. In our code, the Newton-Raphson method was used, provided by the Roots package. 

However, in the code not only the specific cases are developed. We also incorporated the universal formulation of Kepler's equation to solve for the universal anomaly $\chi$ which applied an Stumpff function and the numerical solver.


| **Function Name** | **Inputs** | **Outputs** | **Functionality** |
|-------------------|------------|-------------|-------------------|
| timeSincePeriapsis | stateVector::MyStateVector, mu::Float64, theta::Float64 | time::Float64 | $\theta$ [rad] $\rightarrow \Delta t$ [s] |
| timeSinceTrueAnomaly | stateVector::MyStateVector, mu::Float64, theta1::Float64, theta2::Float64 | time::Float64 | $\Delta t$ between $\theta_1$ & $\theta_2$ |
| stumpffFunction | z::Float64 | S::Float64 | Stumpff function |
| universalAnomaly | stateVector::MyStateVector, mu::Float64, time::Float64 | X::Float64 | $\Delta t \rightarrow \chi$ |
| universalAnomaly\_to\_trueAnomaly | stateVector::MyStateVector, mu::Float64, X::Float64 | theta::Float64 | $\chi \rightarrow \theta$ |
| trueAnomaly\_to\_universalAnomaly | stateVector::MyStateVector, mu::Float64, theta::Float64 | X::Float64 | $\theta \rightarrow \chi$ |
| timeSincePeriapsis\_to\_trueAnomaly | stateVector::MyStateVector, mu::Float64, time::Float64 | theta::Float64 | $\Delta t \rightarrow \theta$ |
| trueAnomaly\_to\_eccentricAnomaly | stateVector::MyStateVector, mu::Float64, theta::Float64 | E::Float64 | $\theta \rightarrow E$ |
| eccentricAnomaly\_to\_trueAnomaly | stateVector::MyStateVector, mu::Float64, E::Float64 | theta::Float64 | $E \rightarrow \theta$ |
| eccentricAnomaly\_to\_meanAnomaly | stateVector::MyStateVector, mu::Float64, E::Float64 | Me::Float64 | $E \rightarrow M_e$ |
| meanAnomaly\_to\_eccentricAnomaly | stateVector::MyStateVector, mu::Float64, Me::Float64 | E::Float64 | $M_e \rightarrow E$ |
| trueAnomaly\_to\_meanAnomaly\_E | stateVector::MyStateVector, mu::Float64, theta::Float64 | Me::Float64 | $\theta \rightarrow M_e$ |
| meanAnomaly\_to\_trueAnomaly\_E | stateVector::MyStateVector, mu::Float64, Me::Float64 | theta::Float64 | $M_e \rightarrow \theta$ |
| universalAnomaly\_to\_eccentricAnomaly | stateVector::MyStateVector, mu::Float64, X::Float64 | E::Float64 | $\chi \rightarrow E$ |
| eccentricAnomaly\_to\_universalAnomaly | stateVector::MyStateVector, mu::Float64, E::Float64 | X::Float64 | $E \rightarrow \chi$ |
| timeSincePeriapsis\_to\_Me | stateVector::MyStateVector, mu::Float64, time::Float64 | Me::Float64 | $\Delta t \rightarrow M_e$ |
| Me\_to\_timeSincePeriapsis | stateVector::MyStateVector, mu::Float64, Me::Float64 | time::Float64 | $M_e \rightarrow \Delta t$ |
| trueAnomaly\_to\_parabolicAnomaly | stateVector::MyStateVector, mu::Float64, theta::Float64 | P::Float64 | $\theta \rightarrow P$ |
| parabolicAnomaly\_to\_trueAnomaly | stateVector::MyStateVector, mu::Float64, P::Float64 | theta::Float64 | $P \rightarrow \theta$ |
| parabolicAnomaly\_to\_meanAnomaly | stateVector::MyStateVector, mu::Float64, P::Float64 | Mp::Float64 | $P \rightarrow M_p$ |
| meanAnomaly\_to\_parabolicAnomaly | stateVector::MyStateVector, mu::Float64, Mp::Float64 | P::Float64 | $M_p \rightarrow P$ |
| trueAnomaly\_to\_meanAnomaly\_P | stateVector::MyStateVector, mu::Float64, theta::Float64 | Mp::Float64 | $\theta \rightarrow M_p$ |
| meanAnomaly\_to\_trueAnomaly\_P | stateVector::MyStateVector, mu::Float64, Mp::Float64 | theta::Float64 | $M_p \rightarrow \theta$ |
| universalAnomaly\_to\_hyperbolicAnomaly | stateVector::MyStateVector, mu::Float64, X::Float64 | H::Float64 | $\chi \rightarrow H$ |
| hyperbolicAnomaly\_to\_universalAnomaly | stateVector::MyStateVector, mu::Float64, H::Float64 | X::Float64 | $H \rightarrow \chi$ |
| timeSincePeriapsis\_to\_Mh | stateVector::MyStateVector, mu::Float64, time::Float64 | Mh::Float64 | $\Delta t \rightarrow M_h$ |
| Mh\_to\_timeSincePeriapsis | stateVector::MyStateVector, mu::Float64, Mh::Float64 | time::Float64 | $M_h \rightarrow \Delta t$ |
| trueAnomaly\_to\_meanAnomaly\_H | stateVector::MyStateVector, mu::Float64, theta::Float64 | Mh::Float64 | $\theta \rightarrow M_h$ |
| meanAnomaly\_to\_trueAnomaly\_H | stateVector::MyStateVector, mu::Float64, Mh::Float64 | theta::Float64 | $M_h \rightarrow \theta$ |
| universalAnomaly\_to\_hyperbolicAnomaly | stateVector::MyStateVector, mu::Float64, X::Float64 | H::Float64 | $\chi \rightarrow H$ |
| hyperbolicAnomaly\_to\_universalAnomaly | stateVector::MyStateVector, mu::Float64, H::Float64 | X::Float64 | $H \rightarrow \chi$ |
| timeSincePeriapsis\_to\_Mh | stateVector::MyStateVector, mu::Float64, time::Float64 | Mh::Float64 | $\Delta t \rightarrow M_h$ |
| Mh\_to\_timeSincePeriapsis | stateVector::MyStateVector, mu::Float64, Mh::Float64 | time::Float64 | $M_h \rightarrow \Delta t$ |


To better understand all the possible directions of conversion between the different anomalies and time, the diagram below can be consulted.

![AnomaliesDiagram](https://user-images.githubusercontent.com/115453770/252145525-2ed4dd46-bdc5-47df-bf61-92051b2e7d7a.png)

Inside the file, the initial orbit determination is also tackled through the Lambert's problem. The functions written in this third table facilitate the determination of the time of flight, conic section parameters, and the solution of Lambert's problem itself.

| **Function Name** | **Inputs** | **Outputs** | **Functionality** |
|-------------------|------------|-------------|-------------------|
| Lambert\_solve | r1::MyVector, r2::MyVector, tF::Float64, mu::Float64, k::MyVector | coe1::MyCOE, coe2::MyCOE | Given two positions $\mathbf{r_1}$ and $\mathbf{r_2}$ and the time of flight between them, find the orbit that connects them |
| Lambert\_conic | r1::MyVector, r2::MyVector, eT::Float64, k::MyVector | coe1::MyCOE, coe2::MyCOE | Obtain the conic section for a given transverse eccentricity $e_T$ using the Avanzini algorithm |
| timeOfFlight | coe1::MyCOE, coe2::MyCOE, mu::Float64, theta2::Float64, mu::Float64 | delta\_t::Float64 | Time of flight between two true anomalies ($\theta_1$ and $\theta_2$) |
