module IdealTwoBodyProblem

using LinearAlgebra
using LinearAlgebraTypes
using Roots

export  MyStateVector, angularMomentumVector, mechanicalEnergy, eccentricityVector, semilatusRectum, periodOrbit, perifocalBasis, periapsis, apoapsis, 
        trajectoryEquation, MyCOE, inclination, semiMajorAxis, orbitType, trueAnomaly, stateVector_to_COE, COE_to_stateVector, 
        escapeVelocity, circularVelocity, velocityPeriapsis, velocityApoapsis, velocityHyperbolicAsymptote, speed_visViva, 
        rightAscensionOfTheAscendingNode, argumentOfPeriapsis, vector_in_ECI_to_perifocal, vector_in_perifocal_to_ECI, 
        timeSincePeriapsis, timeSinceTrueAnomaly, stumpffFunction, universalAnomaly, universalAnomaly_to_trueAnomaly, trueAnomaly_to_universalAnomaly, 
        trueAnomaly_to_eccentricAnomaly, eccentricAnomaly_to_trueAnomaly, eccentricAnomaly_to_meanAnomaly, meanAnomaly_to_eccentricAnomaly,
        trueAnomaly_to_meanAnomaly_E, meanAnomaly_to_trueAnomaly_E, universalAnomaly_to_eccentricAnomaly, eccentricAnomaly_to_universalAnomaly,
        trueAnomaly_to_parabolicAnomaly, parabolicAnomaly_to_trueAnomaly, parabolicAnomaly_to_meanAnomaly,meanAnomaly_to_parabolicAnomaly, 
        trueAnomaly_to_meanAnomaly_P, meanAnomaly_to_trueAnomaly_P, universalAnomaly_to_parabolicAnomaly, parabolicAnomaly_to_universalAnomaly,
        trueAnomaly_to_hyperbolicAnomaly, hyperbolicAnomaly_to_trueAnomaly, hyperbolicAnomaly_to_meanAnomaly, meanAnomaly_to_hyperbolicAnomaly, 
        trueAnomaly_to_meanAnomaly_H, meanAnomaly_to_trueAnomaly_H, universalAnomaly_to_hyperbolicAnomaly, hyperbolicAnomaly_to_universalAnomaly,
        timeSincePeriapsis_to_trueAnomaly, timeSincePeriapsis_to_Me, timeSincePeriapsis_to_Mp, timeSincePeriapsis_to_Mh, Me_to_timeSincePeriapsis,
        Mp_to_timeSincePeriapsis, Mh_to_timeSincePeriapsis, propagate_StateVector, timeOfFlight, Lambert_conic, Lambert_solve

# MyStateVector = Point.position and Point.velocity
struct MyStateVector
    r :: MyVector
    v :: MyVector
end

# Calculate the specfic angular momentum per unit mass of a Point (h = r x v)
function angularMomentumVector(stateVector::MyStateVector)
    # Input 
        # stateVector:: state vector object
    # Output
        # h: mass-specific angular momentum vector
    r = stateVector.r
    v = stateVector.v
    h = cross(r, v)
    return h
end

# Calculate the mass-specific mechanical energy of a point in a reference frame (xi = 0.5*v^2 - mu/r)
function mechanicalEnergy(stateVector::MyStateVector, mu::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
    # Output
        # xi: mass-specific mechanical energy
    r = stateVector.r
    v = stateVector.v
    r_norm = norm(r)
    v_norm = norm(v)
    xi = 0.5*v_norm^2 - mu/r_norm
    return xi
end

# Calculate the eccentricity vector of a point in a reference frame (e = (1/mu)*((v^2 - mu/r)*r - dot(r, v)*v))
function eccentricityVector(stateVector::MyStateVector, mu::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
    # Output
        # e: eccentricity vector
    r = stateVector.r
    v = stateVector.v
    v_norm = norm(v)
    r_norm = norm(r)
    e = (1/mu)*((v_norm^2 - mu/r_norm)*r - dot(r,v)*v)
    return e
end

# Calculate the semiMajorAxis of an orbit (a = -mu/(2*xi))
function semiMajorAxis(stateVector::MyStateVector, mu::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
    # Output
        # a: semi-major axis
    r = stateVector.r
    v = stateVector.v
    xi = mechanicalEnergy(stateVector,mu)
    a = -mu/(2*xi)
    return a
end

# Calculate the semi-latus rectum of an orbit (p = h^2/mu)
function semilatusRectum(stateVector::MyStateVector,mu::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
    # Output
        # p: semilatus rectum
    h = angularMomentumVector(stateVector)
    p = norm(h)^2/mu
    return p
end

# Calculate the period of an orbit (T = 2*pi*sqrt(a^3/mu))
function periodOrbit(stateVector::MyStateVector,mu::Float64,OrbitType::String)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # OrbitType: "Elliptical", "Parabolic", or "Hyperbolic", or "Circular"
    # Output
        # period: period of orbit 
    if OrbitType == "Elliptical"
        a = semiMajorAxis(stateVector,mu)
        period = 2*pi*sqrt(a^3/mu)
    elseif OrbitType == "Parabolic"
        period = Inf
    elseif OrbitType == "Hyperbolic"
        a = semiMajorAxis(stateVector,mu)
        period = 2*pi*sqrt(-a^3/mu)
    elseif OrbitType == "Circular"
        period = 2*pi*sqrt(norm(stateVector.r)^3/mu)
    else
        throw(ArgumentError("OrbitType must be 'Elliptical', 'Parabolic', 'Hyperbolic', or 'Circular'."))
    end       
    return period
end

# Calculate the perifocal basis (p,q,w) from the state vector object and express it in MyBasis
function perifocalBasis(stateVector::MyStateVector, mu::Float64)
    # Input 
        # stateVector: state vector object
        # mu: gravitational parameter
    # Output
        # perifocalBasis: perifocal basis
    r = stateVector.r
    h = angularMomentumVector(stateVector)
    e = eccentricityVector(stateVector,mu)
    if norm(e) == 0
        i = r/norm(r)
    else
        i = e/norm(e)
    end
    k = h/norm(h)
    j = cross(k,i)
    p = i
    q = j
    w = k
    omega = MyVector([0.0,0.0,0.0])
    alpha = MyVector([0.0,0.0,0.0])
    perifocalBasis = MyBasis(p,q,w,omega,alpha)
    return perifocalBasis
end

# Calculate the periapsis 
function periapsis(stateVector::MyStateVector,mu::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
    # Output
        # r_p: periapsis
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    r_p = a*(1-norm(e))
    return r_p
end

# Calculate the apoapsis
function apoapsis(stateVector::MyStateVector,mu::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
    # Output
        # r_a: apoapsis
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    r_a = a*(1+norm(e))
    return r_a
end

# Calculate the escape velocity (v_esc = sqrt(2*mu/r))
function escapeVelocity(r::MyVector,mu::Float64)
    # Input 
        # r: position vector
        # mu: gravitational parameter
    # Output
        # v_esc: escape velocity
    v_esc = sqrt(2*mu/norm(r))
    return v_esc
end

# Calculate the circular velocity (v_circ = sqrt(mu/r))
function circularVelocity(r::MyVector,mu::Float64)
    # Input 
        # r: position vector
        # mu: gravitational parameter
    # Output
        # v_circ: circular velocity
    v_circ = sqrt(mu/norm(r))
    return v_circ
end

function circularVelocity(r_norm::Float64, mu::Float64)
    # Input 
        # r_norm: norm of position vector
        # mu: gravitational parameter
    # Output
        # v_circ: circular velocity
    v_circ = sqrt(mu/r_norm)
    return v_circ
end

# Calculate the velocity at periapsis (v_p = sqrt(mu/a*(1+e)/(1-e)))
function velocityPeriapsis(stateVector::MyStateVector,mu::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
    # Output
        # v_p: velocity at periapsis
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    v_p = sqrt(mu/a*(1+norm(e))/(1-norm(e)))
    return v_p
end

# Calculate the velocity at apoapsis (v_a = sqrt(mu/a*(1-e)/(1+e)))
function velocityApoapsis(stateVector::MyStateVector,mu::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
    # Output
        # v_a: velocity at apoapsis
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    v_a = sqrt(mu/a*(1-norm(e))/(1+norm(e)))
    return v_a
end

# Calculate the velocity at hyperbolic asymptote (v_h = sqrt(-mu/a))
function velocityHyperbolicAsymptote(stateVector::MyStateVector,mu::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
    # Output
        # v_h: velocity at hyperbolic asymptote
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    v_h = sqrt(-mu/a)
    return v_h
end

# Calculate the speed with the Vis-viva equation (v = sqrt(mu*((2/norm(r)) - (1/a))))
function speed_visViva(r_norm::Float64,a::Float64,mu::Float64)
    # Input 
        # r: orbit radius
        # a: semi-major axis
        # mu: gravitational parameter
    # Output
        # v: velocity
    v = sqrt(mu*((2/r_norm) - (1/a)))
    return v
end

# Get the orbit type (Elliptical, Parabolic, Hyperbolic, Circular)
function orbitType(stateVector::MyStateVector,mu::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
    # Output
        # orbitType: orbit type
    e = eccentricityVector(stateVector,mu)
    xi = mechanicalEnergy(stateVector,mu)
    if norm(e) < 1 && xi < 0
        orbitType = "Elliptical"
    elseif norm(e) == 1 && xi == 0
        orbitType = "Parabolic"
    elseif norm(e) > 1 && xi > 0
        orbitType = "Hyperbolic"
    elseif norm(e) == 0 && xi < 0
        orbitType = "Circular"
    else
        throw(ArgumentError("Something went wrong."))
    end
    return orbitType
end

# Calculate the trajectory equation (r = h^2/mu/(1+e*cos(theta))
function trajectoryEquation(stateVector::MyStateVector,mu::Float64,theta::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta: true anomaly
    # Output
        # trajectoryEquation: trajectory equation
    e = eccentricityVector(stateVector,mu)
    h = angularMomentumVector(stateVector)
    trajectoryEquation = norm(h)^2/(mu*(1+norm(e)*cos(theta)))
    return trajectoryEquation
end

# Calculate the inclination angle 
function inclination(stateVector::MyStateVector)
    # Input 
        # stateVector:: state vector object
    # Output
        # inc: inclination angle
    h = angularMomentumVector(stateVector)
    inc = acos(h[3]/norm(h))
    return inc
end

# Calculate the true anomaly (theta)
function trueAnomaly(stateVector::MyStateVector,mu::Float64)
    # Input 
        # stateVector:: state vector object
    # Output
        # theta: true anomaly
    r = stateVector.r 
    perifocal_basis = perifocalBasis(stateVector,mu)
    p = perifocal_basis.i
    q = perifocal_basis.j
    theta = atan(dot(r,q), dot(r,p))
    theta < 0 ? theta = theta + 2*pi : theta = theta

    return theta
end

# Calculate the right ascension of the ascending node (RAAN)
function rightAscensionOfTheAscendingNode(stateVector::MyStateVector)
    # Input 
        # stateVector:: state vector object
    # Output
        # RAAN: right ascension of the ascending node
    h = angularMomentumVector(stateVector)
    k_ECI = [0,0,1]
    node_line = cross(k_ECI, h.c0)
    RAAN = acos(node_line[1]/norm(node_line))
    # if node_line[2] is positive, then RAAN is between 0 and 180 degrees
    # if node_line[2] is negative, then RAAN is between 180 and 360 degrees
    node_line[2] >= 0 ? RAAN = RAAN : RAAN = 2*pi - RAAN

    #RAAN = atan(node_line[2], node_line[1])
    #RAAN < 0 ? RAAN = RAAN + 2*pi : RAAN = RAAN

    return RAAN
end

# Calculate the argument of periapsis (omega)
function argumentOfPeriapsis(stateVector::MyStateVector, mu::Float64)
    # Input 
        # stateVector:: state vector object
    # Output
        # omega: argument of periapsis
    h = angularMomentumVector(stateVector)
    k_ECI = [0,0,1]
    node_line = cross(k_ECI, h.c0)
    e = eccentricityVector(stateVector,mu)
    omega = acos(dot(node_line, e.c0)/(norm(node_line)*norm(e)))

    # if e[3] is positive, then omega is between 0 and 180 degrees
    # if e[3] is negative, then omega is between 180 and 360 degrees
    e[3] >= 0 ? omega = omega : omega = 2*pi - omega

    return omega
end

# MyCOE is a struct that stores all the classical orbital elements
struct MyCOE
    # 3 Euler angles:
    RAAN::Float64           # right ascension of the ascending node (0,2pi)
    inc::Float64              # inclination (0,pi)
    omega::Float64          # argument of periapsis (0,2pi)
    # 2 parameters that give the size and shape of the orbit: (a,e) or (p,e)
    a::Float64              # a: semi-major axis
    p::Float64              # p: semi-latus rectum
    e::Float64              # e: eccentricity
    # 1 parameter that gives the current position of the body in the orbit:
    theta::Float64          # theta: true anomaly (0,2pi)

    function MyCOE(; RAAN::Float64=0.0, inc::Float64=0.0, omega::Float64=0.0, a::Float64=0.0, e::Float64=0.0, theta::Float64=0.0, p::Float64=0.0)
        if a == 0.0 && p != 0.0
            a = p / (1 - e^2)
        elseif p == 0.0 && a != 0.0
            p = a * (1 - e^2)
        elseif p == 0.0 && a == 0.0
            throw(ArgumentError("Either 'p' or 'a' must be provided."))
        else
            throw(ArgumentError("Either 'p' or 'a' must be provided, but not both."))
        end
        new(RAAN, inc, omega, a, p, e, theta)
    end
end

# Calculate the COE (Classical Orbital Elements) from a state vector (r,v)
function stateVector_to_COE(stateVector::MyStateVector,mu::Float64)
    # Input 
        # stateVector: state vector object
        # mu: gravitational parameter
    # Output
        # 3 Euler angles:
            # RAAN: right ascension of the ascending node (0,2pi)
            # inc: inclination (0,pi)
            # omega: argument of periapsis (0,2pi)
        # 2 parameters that give the size and shape of the orbit:
            # a: semi-major axis
            # e: eccentricity
        # 1 parameter that gives the current position of the body in the orbit:
            # theta: true anomaly (0,2pi)
    e = eccentricityVector(stateVector,mu)
    h = angularMomentumVector(stateVector)
    a = semiMajorAxis(stateVector,mu)
    inc = inclination(stateVector)
    if inc == 0 || inc == pi
        # Handle edge cases where inclination is 0 or 180 degrees
        node_line = [1.0, 0.0, 0.0]
        RAAN = 0.0
        omega = atan(e[2], e[1])
    else
        RAAN = rightAscensionOfTheAscendingNode(stateVector)
        omega = argumentOfPeriapsis(stateVector, mu)
    end
    theta = trueAnomaly(stateVector,mu)
    if isapprox(norm(e), 1.0)
        # Handle edge case where eccentricity is 1
        p = semilatusRectum(stateVector,mu)
        coe = MyCOE(RAAN=RAAN, inc=inc, omega=omega, p=p, e=norm(e), theta=theta)
    else
        coe = MyCOE(RAAN=RAAN, inc=inc, omega=omega, a=a, e=norm(e), theta=theta)
    end
    return coe
end

# Calculate the state vector (r,v) from the COE (Classical Orbital Elements)
function COE_to_stateVector(coe::MyCOE,mu::Float64)
    # Input 
        # coe: Classical Orbital Elements
        # mu: gravitational parameter
    # Output
        # r: position vector
        # v: velocity vector        
    # To go from COE to state vector, first compute r and v in the perifocal reference frame vector
    # basis. Then, transform the matrix to the original reference frame in which the state vector
    #should be provided (e.g. ECI-equatorial), by undoing the three rotations associated to the
    # Euler angles.
    h = sqrt(mu*coe.a*(1.0-coe.e^2))
    # r and v in the perifocal reference frame vector basis
    r_PF = h^2/(mu*(1.0+coe.e*cos(coe.theta)))*[cos(coe.theta), sin(coe.theta), 0.0]
    v_PF = mu/h*[-sin(coe.theta), coe.e+cos(coe.theta), 0.0]
    # Transformation from the perifocal reference frame vector basis to the ECI vector basis
    r_ECI = vector_in_perifocal_to_ECI(r_PF, coe)
    v_ECI = vector_in_perifocal_to_ECI(v_PF, coe)
    return MyStateVector(r_ECI,v_ECI)
end 

# Transform a vector from the perifocal reference frame to the ECI reference frame
function vector_in_perifocal_to_ECI(perifocalVector::Vector, coe::MyCOE)
    # Input 
        # perifocalVector: vector in the perifocal reference frame
        # coe: Classical Orbital Elements
    # Output
        # ECI_vector: vector in the ECI reference frame
    
    @assert length(perifocalVector) == 3 "perifocalVector must be a 3-element vector"
    # Transformation matrix from the perifocal reference frame vector basis ECI reference frame
    R_3_omega = [cos(coe.omega) sin(coe.omega) 0.0;-sin(coe.omega) cos(coe.omega) 0.0; 0.0 0.0 1.0]
    R_1_inc = [1.0 0.0 0.0; 0.0 cos(coe.inc) sin(coe.inc); 0.0 -sin(coe.inc) cos(coe.inc)]
    R_3_RAAN = [cos(coe.RAAN) sin(coe.RAAN) 0.0;-sin(coe.RAAN) cos(coe.RAAN) 0.0; 0.0 0.0 1.0]
    T = R_3_RAAN'*R_1_inc'*R_3_omega'
    # vector in the ECI reference frame
    ECI_vector = MyVector(T*perifocalVector)
    return ECI_vector
end

# Transform a vector from the ECI reference frame to the perifocal reference frame
function vector_in_ECI_to_perifocal(ECI_vector::Vector, coe::MyCOE)
    # Input 
        # ECI_vector: vector in the ECI reference frame
        # coe: Classical Orbital Elements
    # Output
        # perifocalVector: vector in the perifocal reference frame
    
    @assert length(ECI_vector) == 3 "ECI_vector must be a 3-element vector"
    # Transformation matrix from the ECI reference frame to the perifocal reference frame
    R_3_RAAN = [cos(coe.RAAN) sin(coe.RAAN) 0.0;-sin(coe.RAAN) cos(coe.RAAN) 0.0; 0.0 0.0 1.0]
    R_1_inc = [1.0 0.0 0.0; 0 cos(coe.inc) sin(coe.inc); 0.0 -sin(coe.inc) cos(coe.inc)]
    R_3_omega = [cos(coe.omega) sin(coe.omega) 0.0;-sin(coe.omega) cos(coe.omega) 0.0; 0.0 0.0 1.0]
    T = R_3_omega*R_1_inc*R_3_RAAN
    # vector in the perifocal reference frame
    perifocalVector = MyVector(T*ECI_vector)
    return perifocalVector
end

# Convert between Anomalies

# Convert from the true anomaly (theta) to the eccentric anomaly (E)
function trueAnomaly_to_eccentricAnomaly(stateVector::MyStateVector,mu::Float64,theta::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta: true anomaly
    # Output
        # E: eccentric anomaly
    e = eccentricityVector(stateVector,mu)
    sinE = sqrt(1-norm(e)^2)*sin(theta)/(1+norm(e)*cos(theta))
    cosE = (norm(e)+cos(theta))/(1+norm(e)*cos(theta))
    E = atan(sinE,cosE)
    return E
end

# Convert from the eccentric anomaly (E) to the true anomaly (theta)
function eccentricAnomaly_to_trueAnomaly(stateVector::MyStateVector,mu::Float64,E::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # E: eccentric anomaly
    # Output
        # theta: true anomaly
    e = eccentricityVector(stateVector,mu)
    sintheta = sqrt(1-norm(e)^2)*sin(E)/(1-norm(e)*cos(E))
    costheta = (cos(E)-norm(e))/(1-norm(e)*cos(E))
    theta = atan(sintheta,costheta)
    return theta
end

# Convert from the eccentric anomaly (E) to the mean anomaly in the elliptic case (Me)
function eccentricAnomaly_to_meanAnomaly(stateVector::MyStateVector,mu::Float64,E::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # E: eccentric anomaly
    # Output
        # Me: mean anomaly
    e = eccentricityVector(stateVector,mu)
    Me = E - norm(e)*sin(E)
    return Me
end

# Convert from the mean anomaly in the elliptic case (Me) to the eccentric anomaly (E)
function meanAnomaly_to_eccentricAnomaly(stateVector::MyStateVector,mu::Float64,Me::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # Me: mean anomaly
    # Output
        # E: eccentric anomaly
    e = eccentricityVector(stateVector,mu)
    f(E) = E - norm(e)*sin(E) - Me
    E = fzero(f, 0.0)
    return E
end

# Convert from the mean anomaly in the elliptic case (Me) to the true anomaly (theta)
function meanAnomaly_to_trueAnomaly_E(stateVector::MyStateVector,mu::Float64,Me::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # Me: mean anomaly
    # Output
        # theta: true anomaly
    E = meanAnomaly_to_eccentricAnomaly(stateVector,mu,Me)
    theta = eccentricAnomaly_to_trueAnomaly(stateVector,mu,E)
    return theta
end

# Convert from true anomaly (theta) to the mean anomaly in the elliptic case (Me)
function trueAnomaly_to_meanAnomaly_E(stateVector::MyStateVector,mu::Float64,theta::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta: true anomaly
    # Output
        # Me: mean anomaly
    E = trueAnomaly_to_eccentricAnomaly(stateVector,mu,theta)
    Me = eccentricAnomaly_to_meanAnomaly(stateVector,mu,E)
    return Me
end

# Convert from true anomaly to the parabolic anomaly (P)
function trueAnomaly_to_parabolicAnomaly(stateVector::MyStateVector,mu::Float64,theta::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta: true anomaly
    # Output
        # P: parabolic anomaly
    P = tan(theta/2)
    return P
end

# Convert from the parabolic anomaly (P) to the true anomaly (theta)
function parabolicAnomaly_to_trueAnomaly(stateVector::MyStateVector,mu::Float64,P::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # P: parabolic anomaly
    # Output
        # theta: true anomaly
    theta = 2*atan(P)
    return theta
end

# Convert from the parabolic anomaly (P) to the mean anomaly in the parabolic case (Mp)
function parabolicAnomaly_to_meanAnomaly(stateVector::MyStateVector,mu::Float64,P::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # P: parabolic anomaly
    # Output
        # Mp: mean anomaly
    Mp = 1/2*P - 1/6*P^3
    return Mp
end

# Convert from the mean anomaly in the parabolic case (Mp) to the parabolic anomaly (P)
function meanAnomaly_to_parabolicAnomaly(stateVector::MyStateVector,mu::Float64,Mp::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # Mp: mean anomaly
    # Output
        # P: parabolic anomaly
    f(P) = 1/2*P - 1/6*P^3 - Mp
    P = fzero(f, 0.0)
    return P
end

# Convert from the mean anomaly in the parabolic case (Mp) to the true anomaly (theta)
function meanAnomaly_to_trueAnomaly_P(stateVector::MyStateVector,mu::Float64,Mp::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # Mp: mean anomaly
    # Output
        # theta: true anomaly
    P = meanAnomaly_to_parabolicAnomaly(stateVector,mu,Mp)
    theta = parabolicAnomaly_to_trueAnomaly(stateVector,mu,P)
    return theta
end

# Convert from true anomaly (theta) to the mean anomaly in the parabolic case (Mp)
function trueAnomaly_to_meanAnomaly_P(stateVector::MyStateVector,mu::Float64,theta::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta: true anomaly
    # Output
        # Mp: mean anomaly
    Mp = 1/2*tan(theta/2)+1/6*tan(theta/2)^3
    return Mp
end

# Convert from true anomaly (theta) to the hyperbolic anomaly (H)
function trueAnomaly_to_hyperbolicAnomaly(stateVector::MyStateVector,mu::Float64,theta::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta: true anomaly
    # Output
        # H: hyperbolic anomaly
    e = eccentricityVector(stateVector,mu)
    H = acosh((norm(e)+cos(theta))/(1+norm(e)*cos(theta)))
    return H
end

# Convert from the hyperbolic anomaly (H) to the true anomaly (theta)
function hyperbolicAnomaly_to_trueAnomaly(stateVector::MyStateVector,mu::Float64,H::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # H: hyperbolic anomaly
    # Output
        # theta: true anomaly
    e = eccentricityVector(stateVector,mu)
    costheta = (cosh(H)-norm(e))/(1-norm(e)*cosh(H))
    sintheta = sqrt(1-norm(e)^2)*sinh(H)/(1-norm(e)*cosh(H))
    theta = atan(sintheta,costheta)
    return theta
end

# Convert from the hyperbolic anomaly (H) to the mean anomaly in the hyperbolic case (Mh)
function hyperbolicAnomaly_to_meanAnomaly(stateVector::MyStateVector,mu::Float64,H::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # H: hyperbolic anomaly
    # Output
        # Mh: mean anomaly
    e = eccentricityVector(stateVector,mu)
    Mh = norm(e)*sinh(H) - H
    return Mh
end

# Convert from the mean anomaly in the hyperbolic case (Mh) to the hyperbolic anomaly (H)
function meanAnomaly_to_hyperbolicAnomaly(stateVector::MyStateVector,mu::Float64,Mh::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # Mh: mean anomaly
    # Output
        # H: hyperbolic anomaly
    f(H) = norm(e)*sinh(H) - H - Mh
    H = fzero(f, 0.0)
    return H
end

# Convert from the mean anomaly in the hyperbolic case (Mh) to the true anomaly (theta)
function meanAnomaly_to_trueAnomaly_H(stateVector::MyStateVector,mu::Float64,Mh::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # Mh: mean anomaly
    # Output
        # theta: true anomaly
    H = meanAnomaly_to_hyperbolicAnomaly(stateVector,mu,Mh)
    theta = hyperbolicAnomaly_to_trueAnomaly(stateVector,mu,H)
    return theta
end

# Convert from true anomaly (theta) to the mean anomaly in the hyperbolic case (Mh)
function trueAnomaly_to_meanAnomaly_H(stateVector::MyStateVector,mu::Float64,theta::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta: true anomaly
    # Output
        # Mh: mean anomaly
    H = trueAnomaly_to_hyperbolicAnomaly(stateVector,mu,theta)
    Mh = hyperbolicAnomaly_to_meanAnomaly(stateVector,mu,H)
    return Mh
end

#=
# Given True Anomaly, find Time Since Periapsis
function timeSincePeriapsis(stateVector::MyStateVector,mu::Float64,theta::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta: true anomaly
    # Output
        # time: time to go from the periapsis to a true anomaly
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    p = semilatusRectum(stateVector,mu)
    if norm(e)==0
        # Circular orbit
        n_c = sqrt(mu/a^3) # mean angular rate
        time = theta/n_c
    elseif norm(e)<1 && norm(e)>0
        # Elliptical orbit
        E = trueAnomaly_to_eccentricAnomaly(stateVector,mu,theta)
        Me = eccentricAnomaly_to_meanAnomaly(stateVector,mu,E)
        n_e = sqrt(mu/a^3) # mean angular rate
        time = Me/n_e 
    elseif norm(e)==1
        # Parabolic orbit
        Mp = trueAnomaly_to_meanAnomaly_P(stateVector,mu,theta)
        n_p = sqrt(mu/(p)^3) # mean angular rate
        time = Mp/n_p
    elseif norm(e)>1
        # Hyperbolic orbit
        H = trueAnomaly_to_hyperbolicAnomaly(stateVector,mu,theta)
        Mh = hyperbolicAnomaly_to_meanAnomaly(stateVector,mu,H)
        n_h = sqrt(-mu/a^3) # mean angular rate
        time = Mh/n_h
    end
    return time
end
=#


# Given Mean Anomaly in the Elliptical Case, find Time Since Periapsis
function Me_to_timeSincePeriapsis(stateVector::MyStateVector,mu::Float64,Me::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # Me: mean anomaly
    # Output
        # time: time to go from the periapsis to a true anomaly
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    n_e = sqrt(mu/a^3) # mean angular rate
    time = Me/n_e 
    return time
end

# Given Mean Anomaly in the Parabolic Case, find Time Since Periapsis
function Mp_to_timeSincePeriapsis(stateVector::MyStateVector,mu::Float64,Mp::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # Mp: mean anomaly
    # Output
        # time: time to go from the periapsis to a true anomaly
    p = semilatusRectum(stateVector,mu)
    n_p = sqrt(mu/(p)^3) # mean angular rate
    time = Mp/n_p
    return time
end

# Given Mean Anomaly in the Hyperbolic Case, find Time Since Periapsis
function Mh_to_timeSincePeriapsis(stateVector::MyStateVector,mu::Float64,Mh::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # Mh: mean anomaly
    # Output
        # time: time to go from the periapsis to a true anomaly
    a = semiMajorAxis(stateVector,mu)
    n_h = sqrt(-mu/a^3) # mean angular rate
    time = Mh/n_h
    return time
end

# Given Time Since Periapsis, find the Mean Anomaly in the Elliptical Case
function timeSincePeriapsis_to_Me(stateVector::MyStateVector,mu::Float64,time::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # time: time to go from the periapsis to a true anomaly
    # Output
        # Me: mean anomaly
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    n_e = sqrt(mu/a^3) # mean angular rate
    Me = n_e*time
    return Me
end 

# Given Time Since Periapsis, find the Mean Anomaly in the Parabolic Case
function timeSincePeriapsis_to_Mp(stateVector::MyStateVector,mu::Float64,time::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # time: time to go from the periapsis to a true anomaly
    # Output
        # Mp: mean anomaly
    p = semilatusRectum(stateVector,mu)
    n_p = sqrt(mu/(p)^3) # mean angular rate
    Mp = n_p*time
    return Mp
end

# Given Time Since Periapsis, find the Mean Anomaly in the Hyperbolic Case
function timeSincePeriapsis_to_Mh(stateVector::MyStateVector,mu::Float64,time::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # time: time to go from the periapsis to a true anomaly
    # Output
        # Mh: mean anomaly
    a = semiMajorAxis(stateVector,mu)
    n_h = sqrt(-mu/a^3) # mean angular rate
    Mh = n_h*time
    return Mh
end


# Given True Anomaly, find Time Since Periapsis (t_0 = 0)
function timeSincePeriapsis(stateVector::MyStateVector,mu::Float64,theta::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta: true anomaly
    # Output
        # time: time to go from the periapsis to a true anomaly
    X = trueAnomaly_to_universalAnomaly(stateVector,mu,theta)
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    time = (norm(e)*X^3*stumpffFunction(X^2/a)+a*(1-norm(e))*X)/sqrt(mu)
    return time
end


# Given 2 True Anomalies, find Time Since the 1st True Anomaly
function timeSinceTrueAnomaly(stateVector::MyStateVector,mu::Float64,theta1::Float64,theta2::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta1: true anomaly 1
        # theta2: true anomaly 2
    # Output
        # time: time to go from theta1 to theta2
    time = timeSincePeriapsis(stateVector,mu,theta2) - timeSincePeriapsis(stateVector,mu,theta1)
    period = periodOrbit(stateVector,mu,orbitType(stateVector,mu))
    time < 0 ? time = time + period : time = time
    return time
end

# Stumpff function (Solve Kepler's equation)
function stumpffFunction(z::Float64)
    # Input 
        # z: Stumpff variable
    # Output
        # S(z): Stumpff equation
    if z > 0
        S = (sqrt(z) - sin(sqrt(z)))/sqrt(z^3)
    elseif z < 0
        S = (sinh(sqrt(-z)) - sqrt(-z))/sqrt((-z)^3)
    else
        S = 1/6
    end
    return S
end

# Caculate the Universal Anomaly X (Solve Kepler's equation) knowing the time since the periapsis, use Stumpff equation
# to solve the nonlinear equation sqrt(mu)*(t-t0) = e*X^3*S(X^2/a)+a*(1-e)*X
function universalAnomaly(stateVector::MyStateVector,mu::Float64,time::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # time: time to go from the periapsis to a true anomaly (time = t-t0)
    # Output
        # X: Universal Anomaly
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    # Solve the nonlinear equation with Roots
    f(X) = sqrt(mu)*(time) - norm(e)*X^3*stumpffFunction(X^2/a) - a*(1-norm(e))*X
    X = fzero(f, 0.0)
    return X
end

# Convert from the Universal Anomaly (X) to the Eccentric Anomaly (E)
function universalAnomaly_to_eccentricAnomaly(stateVector::MyStateVector,mu::Float64,X::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # X: Universal Anomaly
    # Output
        # E: Eccentric Anomaly
    a = semiMajorAxis(stateVector,mu)
    E = X/sqrt(a)
    return E
end

# Convert from the Eccentric Anomaly (E) to the Universal Anomaly (X)
function eccentricAnomaly_to_universalAnomaly(stateVector::MyStateVector,mu::Float64,E::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # E: Eccentric Anomaly
    # Output
        # X: Universal Anomaly
    a = semiMajorAxis(stateVector,mu)
    X = sqrt(a)*E
    return X
end

# Convert from the Universal Anomaly (X) to the Parabolic Anomaly (P)
function universalAnomaly_to_parabolicAnomaly(stateVector::MyStateVector,mu::Float64,X::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # X: Universal Anomaly
    # Output
        # P: Parabolic Anomaly
    p = semilatusRectum(stateVector,mu)
    P = X/sqrt(p)
    return P
end

# Convert from the Parabolic Anomaly (P) to the Universal Anomaly (X)
function parabolicAnomaly_to_universalAnomaly(stateVector::MyStateVector,mu::Float64,P::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # P: Parabolic Anomaly
    # Output
        # X: Universal Anomaly
    p = semilatusRectum(stateVector,mu)
    X = sqrt(p)*P
    return X
end

# Convert from the Universal Anomaly (X) to the Hyperbolic Anomaly (H)
function universalAnomaly_to_hyperbolicAnomaly(stateVector::MyStateVector,mu::Float64,X::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # X: Universal Anomaly
    # Output
        # H: Hyperbolic Anomaly
    a = semiMajorAxis(stateVector,mu)
    H = X/sqrt(-a)
    return H
end

# Convert from the Hyperbolic Anomaly (H) to the Universal Anomaly (X)
function hyperbolicAnomaly_to_universalAnomaly(stateVector::MyStateVector,mu::Float64,H::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # H: Hyperbolic Anomaly
    # Output
        # X: Universal Anomaly
    a = semiMajorAxis(stateVector,mu)
    X = sqrt(-a)*H
    return X
end

# Convert from the Universal Anomaly (X) to the True Anomaly (theta)
function universalAnomaly_to_trueAnomaly(stateVector::MyStateVector,mu::Float64,X::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # X: Universal Anomaly
    # Output
        # theta: True Anomaly
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    if norm(e)== 0 
        # Circular orbit
        theta = X/sqrt(a)
    elseif norm(e) < 1 && norm(e) > 0
        # Elliptical orbit
        E = universalAnomaly_to_eccentricAnomaly(stateVector,mu,X)
        theta = eccentricAnomaly_to_trueAnomaly(stateVector,mu,E)
    elseif norm(e) == 1
        # Parabolic orbit
        P = universalAnomaly_to_parabolicAnomaly(stateVector,mu,X)
        theta = parabolicAnomaly_to_trueAnomaly(stateVector,mu,P)
    elseif norm(e) > 1
        # Hyperbolic orbit
        H = universalAnomaly_to_hyperbolicAnomaly(stateVector,mu,X)
        theta = hyperbolicAnomaly_to_trueAnomaly(stateVector,mu,H)
    end
    return theta
end

# Convert from True Anomaly to Universal Anomaly
function trueAnomaly_to_universalAnomaly(stateVector::MyStateVector,mu::Float64,theta::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # theta: True Anomaly
    # Output
        # X: Universal Anomaly
    e = eccentricityVector(stateVector,mu)
    a = semiMajorAxis(stateVector,mu)
    if norm(e)== 0 
        # Circular orbit
        X = sqrt(a)*theta
    elseif norm(e) < 1 && norm(e) > 0
        # Elliptical orbit
        E = trueAnomaly_to_eccentricAnomaly(stateVector,mu,theta)
        X = eccentricAnomaly_to_universalAnomaly(stateVector,mu,E)
    elseif norm(e) == 1
        # Parabolic orbit
        P = trueAnomaly_to_parabolicAnomaly(stateVector,mu,theta)
        X = parabolicAnomaly_to_universalAnomaly(stateVector,mu,P)
    elseif norm(e) > 1
        # Hyperbolic orbit
        H = trueAnomaly_to_hyperbolicAnomaly(stateVector,mu,theta)
        X = hyperbolicAnomaly_to_universalAnomaly(stateVector,mu,H)
    end
    return X
end

# Given Time Since Periapsis, find True Anomaly
function timeSincePeriapsis_to_trueAnomaly(stateVector::MyStateVector,mu::Float64,time::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # time: time to go from the periapsis to a true anomaly
    # Output
        # theta: True Anomaly
    X = universalAnomaly(stateVector,mu,time)
    theta = universalAnomaly_to_trueAnomaly(stateVector,mu,X)
    return theta
end

# Propagate the state vector (r,v) a desired time (delta-t) into the future
function propagate_StateVector(stateVector::MyStateVector,mu::Float64,delta_t::Float64)
    # Input 
        # stateVector:: state vector object
        # mu: gravitational parameter
        # delta_t: time to propagate the state vector
    # Output
        # stateVector2: propagated state vector
    # Calculate the true anomaly
    coe1 = stateVector_to_COE(stateVector,mu)
    theta1 = coe1.theta
    # Calculate the time since the periapsis
    t1 = timeSincePeriapsis(stateVector,mu,theta1)
    # Calculate the new time since the periapsis
    tf = t1 + delta_t
    # Calculate the new true anomaly
    X = universalAnomaly(stateVector,mu,tf)
    theta2 = universalAnomaly_to_trueAnomaly(stateVector,mu,X)
    # Calculate the new state vector
    if isapprox(coe1.e,1.0)
        coe2 = MyCOE(RAAN=coe1.RAAN,inc=coe1.inc,omega=coe1.omega,p=coe1.p,e=coe1.e,theta=theta2)
    else
        coe2 = MyCOE(RAAN=coe1.RAAN,inc=coe1.inc,omega=coe1.omega,a=coe1.a,e=coe1.e,theta=theta2)
    end
    stateVector2 = COE_to_stateVector(coe2,mu)
    return stateVector2
end

# Calculate the time of flight
function timeOfFlight(coe1::MyCOE,coe2::MyCOE,mu::Float64)
    # Input 
        # coe1: Classical Orbital Elements at the initial time
        # coe2: Classical Orbital Elements at the final time
        # mu: gravitational parameter
    # Output
        # delta_t: time of flight
    # Calculate the true anomaly
    theta1 = coe1.theta
    theta2 = coe2.theta
    # Obtain the state vector
    stateVector1 = COE_to_stateVector(coe1,mu)
    stateVector2 = COE_to_stateVector(coe2,mu)
    # Calculate the time since the periapsis
    t1 = timeSincePeriapsis(stateVector1,mu,theta1)
    t2 = timeSincePeriapsis(stateVector2,mu,theta2)
    # Calculate the time of flight
    if t2 < t1
        delta_t = t2 + periodOrbit(stateVector1,mu, orbitType(stateVector1,mu)) - t1
    else
        delta_t = t2 - t1
    end
    return delta_t
end

# Calculate the time of flight
function timeOfFlight(coe1::MyCOE,theta2::Float64,mu::Float64)
    # Input 
        # coe1: Classical Orbital Elements at the initial time
        # theta2: true anomaly at the final time
        # mu: gravitational parameter
    # Output
        # delta_t: time of flight
    # Calculate the initial true anomaly
    theta1 = coe1.theta
    # Obtain the state vector
    stateVector1 = COE_to_stateVector(coe1,mu)
    # Calculate the time since the periapsis
    t1 = timeSincePeriapsis(stateVector1,mu,theta1)
    t2 = timeSincePeriapsis(stateVector1,mu,theta2)
    # Calculate the time of flight
    if t2 < t1
        delta_t = t2 + periodOrbit(stateVector1,mu, orbitType(stateVector1,mu)) - t1
    else
        delta_t = t2 - t1
    end
    return delta_t
end

# LAMBERT's PROBLEM

# Given two position vectors r1 and r2 of a body and the time of flight between them (tF), find the orbit that connects them
function Lambert_solve(r1::MyVector,r2::MyVector,tF::Float64,mu::Float64,k::MyVector) 
    # Input:
        # r1: initial position vector
        # r2: final position vector
        # tF: time of flight
        # mu: gravitational parameter
        # K: plane normal (only used when rr1 and rr2 are collinear)
    # Output:
       # coe1: Classical Orbital Elements of the connecting orbit at the initial time
       # coe2: Classical Orbital Elements of the connecting orbit at the final time

    tol = 1.0e-8
    maxiter = 1000
    # Calculate the fundamental eccentricity
    r1_norm = norm(r1)
    r2_norm = norm(r2)
    c = norm(r2 - r1)
    eF = -(r2_norm - r1_norm)/c
    eTP = sqrt(1.0 - eF^2)
    
    # Solve for the transverse eccentricity using the newton's method
    f(eT) = error_function(mu,r1,r2,eT,tF,k)
    eT = fzero(f, eTP*0.999, xtol=tol, maxevals=maxiter)


    coe1, coe2 = Lambert_conic(r1,r2,eT,k)
    return coe1, coe2
    
end


function Lambert_conic(r1::MyVector, r2::MyVector, eT::Float64, k::MyVector)
    # Input:
    # r1: initial position vector
    # r2: final position vector
    # eT: transverse eccentricity
    # k: plane normal (only used when r1 and r2 are collinear)
    
    # Output:
    # coe1: Classical Orbital Elements of the first position in the orbit
    # coe2: Classical Orbital Elements of the second position in the orbit
    
    r1_norm = norm(r1)
    r2_norm = norm(r2)
    
    # Calculate the chord vector
    c = r2 - r1
    c_norm = norm(c)
    i = c / c_norm
    
    # Compute the j and k vectors
    temp = cross(r1, r2)
    
    if norm(temp) < 1e-8
        k = k / norm(k)
    else
        # The sign of k*temp forces k vector to point in the direction of the input
        s = sign(dot(k, temp))
        
        if s != 0 # Just for safety
            k = s * temp / norm(temp)
        else
            k = temp / norm(temp)
        end
    end
    
    j = cross(k, i)
    
    # Compute the angle
    Dphi = atan(dot(r2, j), dot(r2, i)) - atan(dot(r1, j), dot(r1, i))
    
    # Compute the parameters of the fundamental ellipse
    eF = -(r2_norm - r1_norm) / c_norm
    aF = (r1_norm + r2_norm) / 2
    pF = aF * (1.0 - eF^2)
    
    # Compute the eccentricity vector and scalar eccentricity
    e_vec = eF * i + eT * j
    e = norm(e_vec)
    
    # Compute the semi-major axis
    if abs(e - 1.0) < 1e-8 # Parabolic orbit
        p = pF
        a = Inf
    elseif e < 1.0 # Elliptic orbit
        p = pF - eT * r1_norm * r2_norm / c_norm * sin(Dphi)
        a = p / (1.0 - e^2)
    else # Hyperbolic orbit
        p = abs(pF - eT * r1_norm * r2_norm / c_norm * sinh(Dphi))
        a = p / (1.0 - e^2)
    end
    
    # Compute the inclination and RAAN
    temp = cross(MyVector([0.0, 0.0, 1.0]), k) # Line of nodes vector
    
    if norm(temp) < 1e-8
        n = MyVector([1.0, 0.0, 0.0])
        m = MyVector([0.0, 1.0, 0.0])
        inc = 0.0
        RAAN = 0.0
    else
        n = temp / norm(temp) # Vector along the line of nodes
        m = cross(k, n) # Vector in the orbital plane (perpendicular to the line of nodes)
        RAAN = mod(atan(n[2], n[1]), 2 * pi)
        inc = acos(k[3])
    end
    
    # Compute the argument of periapsis
    if e < 1e-8 || abs(e - 1.0) < 1e-8
        omega = 0.0
    else
        omega = mod(atan(dot(e_vec, m), dot(e_vec, n)), 2 * pi)
    end
    
    # Compute the true anomalies (this is only valid for elliptic orbits)
    theta1 = mod(atan(dot(r1, m), dot(r1, n)) - omega, 2 * pi)
    theta2 = mod(atan(dot(r2, m), dot(r2, n)) - omega, 2 * pi)

    if e > 1.0 
        # Restrict the range to 0 theta < theta_asymptote
        theta_asymptote = acos(-1.0 / e)
        if theta1 > theta_asymptote
            theta1 = theta_asymptote
        end
        if theta2 > theta_asymptote
            theta2 = theta_asymptote
        end
    end 
    
    # Return the COEs
    if abs(e - 1.0) < 1e-8
        coe1 = MyCOE(RAAN = RAAN, inc = inc, omega = omega, p = p, e = e, theta = theta1)
        coe2 = MyCOE(RAAN = RAAN, inc = inc, omega = omega, p = p, e = e, theta = theta2)
    else
        coe1 = MyCOE(RAAN = RAAN, inc = inc, omega = omega, a = a, e = e, theta = theta1)
        coe2 = MyCOE(RAAN = RAAN, inc = inc, omega = omega, a = a, e = e, theta = theta2)
    end
    
    return coe1, coe2
end

function error_function(mu, r1, r2, eT, tf, k)
    coe1, coe2 = Lambert_conic(r1, r2, eT, k)
    err = tf - timeOfFlight(coe1, coe2, mu)
    return err
end



end # module