# IDEAL TWO BODY PROBLEM TESTS

using LinearAlgebra
using AstrodynamicsEdu.LinearAlgebraTypes
using AStrodynamicsEdu.IdealTwoBodyProblem
using Test

# Test MyStateVector constructor
function test_MyStateVector()
    r1 = MyVector(1.0, 2.0, 3.0)
    v1 = MyVector(1.0, 2.0, 3.0)
    state1 = MyStateVector(r1, v1)
    @test state1.r == r1
    @test state1.v == v1
    # Test MyStateVector object creation with invalid input
    @test_throws MethodError MyStateVector(1, 2)
    @test_throws MethodError MyStateVector(r1,[1.0,2.0,3.0])
    @test_throws MethodError MyStateVector()
end

# Test MyCOE constructor
function test_MyCOE()
    # Test MyCOE object creation with valid input
    RAAN1 = 1.0
    inc1 = 0.5
    omega1 = 1.2
    e1 = 0.5
    theta1 = 6.0
    a1 = 10000.0
    p1  = a1*(1.0 - e1^2)
    # If the semi-major axis is provided
    coe1 = MyCOE(RAAN=RAAN1, inc=inc1, omega=omega1, a=a1, e=e1, theta=theta1)
    @test coe1.RAAN == RAAN1
    @test coe1.inc == inc1
    @test coe1.omega == omega1
    @test coe1.a == a1
    @test coe1.e == e1
    @test coe1.theta == theta1
    @test coe1.p == p1
    # If the semilatus rectum is provided
    coe2 = MyCOE(RAAN=RAAN1, inc=inc1, omega=omega1, p=p1, e=e1, theta=theta1)
    @test coe2.RAAN == RAAN1
    @test coe2.inc == inc1
    @test coe2.omega == omega1
    @test coe2.a == a1
    @test coe2.e == e1
    @test coe2.theta == theta1
    @test coe2.p == coe1.p
    @test coe2.a == coe1.a
    # Even if given in a different order, the constructor should still work
    coe3 = MyCOE(inc=inc1, RAAN=RAAN1, omega=omega1, p=p1, e=e1, theta=theta1)
    @test coe3.RAAN == RAAN1
    @test coe3.inc == inc1
    @test coe3.omega == omega1
    @test coe3.a == a1
    @test coe3.e == e1
    @test coe3.theta == theta1
    @test coe3.p == p1
    # Test MyCOE object creation with invalid input
    @test_throws MethodError MyCOE(1, 2, 3, 4, 5, 6)
    @test_throws ArgumentError MyCOE(RAAN=1.0, inc=0.5, omega=1.2, e=0.5, theta=6.0)
    @test_throws ArgumentError MyCOE(RAAN=1.0, inc=0.5, omega=1.2, e=0.5, theta=6.0, a=10000.0, p=8000.0)
    @test_throws ArgumentError MyCOE()
end


# Test individual functions
function test_I2BPfunctions()
    println("EXAMPLE 1")
    # Data:
    r = MyVector([-5200,-3500, 4637]) # km/s
    v = MyVector([5.2,6.4,2.5]) # km/s
    mu = 3.986e5 # km^3/s^2
    stateVector = MyStateVector(r,v)
    perifocal_basis = perifocalBasis(stateVector,mu)
    println(perifocal_basis)
    period = periodOrbit(stateVector,mu,orbitType(stateVector,mu))
    println("Period: ", period, " s")
    v_perigee = velocityPeriapsis(stateVector,mu)
    println("Velocity at perigee: ", v_perigee, " km/s")
    v_apogee = velocityApoapsis(stateVector,mu)
    println("Velocity at apogee: ", v_apogee, " km/s")
    v_escape = escapeVelocity(r,mu)
    println("Orbit type: ", orbitType(stateVector,mu))
    coe = stateVector_to_COE(stateVector,mu)
    println("Trajectory equation: ", trajectoryEquation(stateVector,mu,coe.theta))

    println("------------------------------------")
    println("EXAMPLE 2")
    # Data:
    a = 10000.0 #km
    e = 0.17
    inc = deg2rad(45.0) #rad
    Omega = deg2rad(120.0) #rad
    omega = deg2rad(15.0) #rad
    theta = deg2rad(120.0) #rad
    mu = 4.0e5 # Earth's gravitational parameter
    coe2 = MyCOE(RAAN=Omega, inc=inc, omega=omega, a=a, e=e, theta=theta)
    # Calculate:
    stateVector2 = COE_to_stateVector(coe2,mu)
    println("Position vector: ", stateVector2.r, "km")
    println("Velocity vector: ", stateVector2.v, "km/s")

    # Problem 2: Kepler’s Equation.
    # A satellite is in an elliptical orbit around the earth. At perigee, the satellite is at 600 km altitude and moving 9.0 km/s.
    # a. Find the following: Semi-major axis. Magnitude of Eccentricity. Distance to apogee. Orbital period
    # b. Using Kepler’s Equation, determine the time, T70, for the satellite to move from perigee to a true anomaly of θ =80°.
    println("------------------------------------")
    println("PROBLEM 2")
    # Data: 
    r = MyVector([6378.0+600.0,0.0,0.0]) # km
    v = MyVector([0.0,9.0,0.0]) # km/s
    stateVector = MyStateVector(r,v)
    mu = 4.00e5 # km^3/s^2
    theta1 = deg2rad(80.0) # rad
    # Calculate:
    coe = stateVector_to_COE(stateVector,mu)
    println("Semi-major axis: ", coe.a, " km")
    println("Magnitude of Eccentricity: ", coe.e)
    println("Distance to apogee: ", coe.a*(1+coe.e), " km")
    period = periodOrbit(stateVector,mu,orbitType(stateVector,mu))
    println("Orbital period: ", period, " s")
    time = timeSincePeriapsis(stateVector,mu,theta1)
    println("Time to go from the periapsis to a true anomaly of 80 degrees: ", time, " s")
    # Now, given the time, calculate the true anomaly:
    theta2 = timeSincePeriapsis_to_trueAnomaly(stateVector,mu,time)
    println("True anomaly: ", theta2*180/pi, " degrees")
    println("----")
    X = trueAnomaly_to_universalAnomaly(stateVector,mu,theta1)
    println("Universal anomaly: ", X)
    thetaFinal = deg2rad(75)
    coeFinal = MyCOE(RAAN=coe.RAAN, inc=coe.inc, omega=coe.omega, a=coe.a, e=coe.e, theta=thetaFinal)
    time2 = timeSincePeriapsis(stateVector,mu,thetaFinal)
    println("Time to go from the periapsis to a true anomaly of 75 degrees: ", time2, " s")
    coe1 = MyCOE(RAAN=coe.RAAN, inc=coe.inc, omega=coe.omega, a=coe.a, e=coe.e, theta=theta1)
    tof = timeOfFlight(coe1,thetaFinal,mu)
    println("Time of flight between 80º and 75º: ", tof, " s")
    tof = timeOfFlight(coe1,coeFinal,mu)
    println("Time of flight between 80º and 75º: ", tof, " s")
    time_theta1_thetaFinal = timeSinceTrueAnomaly(stateVector,mu,theta1,thetaFinal)
    println("Time to go from theta1 to thetaFinal: ", time_theta1_thetaFinal, " s")

end

# Homework 5. Problem 1. Orbit Elements (ELLIPTIC ORBIT)
function test_OrbitElements()
    # Data:
    r = MyVector([-5200,-3500, 4637]) # km/s
    v = MyVector([5.2,6.4,2.5]) # km/s
    mu = 3.986e5 # km^3/s^2
    stateVector = MyStateVector(r,v)
    # Calculate:
    h = angularMomentumVector(stateVector)
    @test isapprox(h,MyVector([-3.8427e4, 3.7112e4, -1.5080e4]), rtol=1e-3, atol=1e-4)
    e = eccentricityVector(stateVector,mu)
    @test isapprox(e,MyVector([0.1919, 0.4046, 0.5067]), rtol=1e-3, atol=1e-4)
    @test isapprox(norm(e),0.6762, rtol=1e-3, atol=1e-4)
    a = semiMajorAxis(stateVector,mu) 
    @test isapprox(a,14237.0,rtol=1e-3, atol=1e-4)
    inc = inclination(stateVector)
    @test isapprox(rad2deg(inc),105.76,rtol=1e-4, atol=1e-4)
    RAAN = rightAscensionOfTheAscendingNode(stateVector)
    @test isapprox(rad2deg(RAAN),226.0,rtol=1e-4, atol=1e-4)
    omega = argumentOfPeriapsis(stateVector,mu)
    @test isapprox(rad2deg(omega),128.9,rtol=1e-4, atol=1e-4) # degrees
    theta = trueAnomaly(stateVector,mu)
    @test isapprox(rad2deg(theta),269.3,rtol=1e-4, atol=1e-4) # degrees
    coe = stateVector_to_COE(stateVector,mu)
    p = semilatusRectum(stateVector,mu)
    @test coe.a == a
    @test coe.e == norm(e)
    @test coe.inc == inc
    @test coe.RAAN == RAAN
    @test coe.omega == omega
    @test coe.theta == theta
    @test isapprox(coe.p, p) 
end

# Homework 5. Problem 2. State Vector (ELLIPTIC ORBIT)
function test_COE2StateVector()
    # Data:
    a = 10000.0 #km
    e = 0.17
    inc = deg2rad(45.0) #rad
    Omega = deg2rad(120.0) #rad
    omega = deg2rad(15.0) #rad
    theta = deg2rad(120.0) #rad
    mu = 3.986e5 # Earth's gravitational parameter
    coe2 = MyCOE(RAAN=Omega, inc=inc, omega=omega, a=a, e=e, theta=theta)
    # Calculate:
    stateVector2 = COE_to_stateVector(coe2,mu)
    @test isapprox(stateVector2.r,MyVector([-843.3, -9153.0, 5307.0]),rtol=1e-3, atol=1e-4) # km
    @test isapprox(stateVector2.v,MyVector([4.536, -2.938, -2.460]),rtol=1e-3, atol=1e-4) # km/s
end

# Homework 3. Problem 2. Kepler’s Equation. (ELLIPTIC ORBIT)
function test_KeplerElliptic()
    # Problem 2: Kepler’s Equation.
    # A satellite is in an elliptical orbit around the earth. At perigee, the satellite is at 600 km altitude and moving 9.0 km/s.
    # Data: 
    rp = 6378.0+600.0 # km
    vp = 9.0 # km/s
    mu = 4.00e5 # km^3/s^2
    r = MyVector([6378.0+600.0,0.0,0.0]) # km
    v = MyVector([0.0,9.0,0.0]) # km/s
    stateVector = MyStateVector(r,v)
    theta1 = deg2rad(80.0) # rad
    # Calculate:
    # a. Find the following: Semi-major axis. Magnitude of Eccentricity. Distance to apogee. Orbital period
    coe = stateVector_to_COE(stateVector,mu)
    @test isapprox(coe.a, 11833.6, rtol=1e-2, atol=1e-4) # km
    @test isapprox(coe.e, 0.411, rtol=1e-2, atol=1e-4)
    rA = coe.a*(1+coe.e) # km
    @test isapprox(rA, 16697.2, rtol=1e-2, atol=1e-4) # km
    period = periodOrbit(stateVector,mu,orbitType(stateVector,mu))
    @test isapprox(period, 12788.7, rtol=1e-2, atol=1e-4) # s
    # b. Using Kepler’s Equation, determine the time, T70, for the satellite to move from perigee to a true anomaly of θ =80°.
    time = timeSincePeriapsis(stateVector,mu,theta1)
    @test isapprox(time, 1321.4, rtol=1e-2, atol=1e-4) # s
    # Now, given the time, calculate the true anomaly:
    theta2 = timeSincePeriapsis_to_trueAnomaly(stateVector,mu,time)
    @test isapprox(theta2,theta1, rtol=1e-3, atol=1e-4) # rad
    # Time to go from theta = 80º to theta = 75º
    thetaFinal = deg2rad(75)
    coeFinal = MyCOE(RAAN=coe.RAAN, inc=coe.inc, omega=coe.omega, a=coe.a, e=coe.e, theta=thetaFinal)
    time2 = timeSincePeriapsis(stateVector,mu,thetaFinal)
    time2thetaFinal = period + time2 - time
    coe1 = MyCOE(RAAN=coe.RAAN, inc=coe.inc, omega=coe.omega, a=coe.a, e=coe.e, theta=theta1)
    tof1 = timeOfFlight(coe1,thetaFinal,mu)
    @test isapprox(time2thetaFinal, tof1, rtol=1e-2, atol=1e-4) # s
    tof2 = timeOfFlight(coe1,coeFinal,mu)
    @test tof2 == tof1
    time_theta1_thetaFinal = timeSinceTrueAnomaly(stateVector,mu,theta1,thetaFinal)
    @test isapprox(time_theta1_thetaFinal, tof, rtol=1e-4, atol=1e-4) # s

    # Check the conversion between anomalies
    E = trueAnomaly_to_eccentricAnomaly(stateVector,mu,theta1)
    @test isapprox(E, 0.9936, rtol=1e-2, atol=1e-4) # rad
    theta2 = eccentricAnomaly_to_trueAnomaly(stateVector,mu,E)
    @test isapprox(theta2,theta1) # rad
    Me = eccentricAnomaly_to_meanAnomaly(stateVector,mu,E)
    @test isapprox(Me, 0.649, rtol=1e-2, atol=1e-4) # rad
    E_2 = meanAnomaly_to_eccentricAnomaly(stateVector,mu,Me)
    @test isapprox(E_2, E) # rad
    theta3 = meanAnomaly_to_trueAnomaly_E(stateVector,mu,Me)
    @test isapprox(theta3,theta1) # rad

end



# Run the tests
@testset "IdealTwoBodyProblem tests" begin 
    test_MyStateVector()
    test_MyCOE()
    #test_I2BPfunctions()
    test_OrbitElements()
    test_COE2StateVector()
    test_KeplerElliptic()
end

