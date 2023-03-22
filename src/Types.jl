module Types

using LinearAlgebra

export VectorBasis, CanonicalVectorBasis, CylindricalVectorBasis, SphericalVectorBasis, IntrinsicVectorBasis, VectorWithBasis, Point, ReferenceFrame, Particle
Particle

# Vector Basis Type defined by 3 orthonormal vectors in a 3x3 matrix, 1 angular velocity vector and 1 angular acceleration vector
mutable struct VectorBasis
    basis::Matrix{Real} # 3x3 orthonormal matrix
    omega::Vector{Real} # angular velocity vector
    alpha::Vector{Real} # angular acceleration vector
end

# Define a constructor that takes in the unitary vectors, the angular velocity vector, and the angular acceleration vector
function VectorBasisConstructor(basis::Matrix{Real}, omega::Vector{Real}, alpha::Vector{Real})
    # Check that the input matrices and vectors have the correct dimensions
    @assert size(basis) == (3, 3) "Unitary vectors must be a 3x3 matrix"
    @assert length(omega) == 3 "Angular velocity vector must have length 3"
    @assert length(alpha) == 3 "Angular acceleration vector must have length 3"
    
    # Check that the unitary vectors are orthonormal
    @assert isapprox(det(basis), 1.0, atol=eps()) "Unitary vectors must be orthonormal"
    
    # Construct the VectorBasis object and return it
    return VectorBasis(basis, omega, alpha)
end

# Define a default constructor that initializes the basis to the Canonical Vector Basis
VectorBasis() = VectorBasis(Matrix{Float64}(I, 3, 3), [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

# Define the Canonical Vector Basis as a VectorBasis object
CanonicalVectorBasis = VectorBasis()

# Cylindrical Vector Basis that is an instance of VectorBasis
theta = 0.0 # Define theta (angle between the x-axis and the radial vector)

CylindricalVectorBasis = VectorBasis(
    Float64[cos(theta) sin(theta) 0.0; -sin(theta) cos(theta) 0.0; 0.0 0.0 1.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0]
)

# Spherical Vector Basis that is an instance of VectorBasis
    theta = 0.0 # Define theta (angle between the x-axis and the radial vector)
    phi = 0.0 # Define phi (angle between the z-axis and the radial vector)

SphericalVectorBasis = VectorBasis(
    Float64[sin(phi)*cos(theta) sin(phi)*sin(theta) cos(phi);
        cos(phi)*cos(theta) cos(phi)*sin(theta) -sin(phi);
        -sin(theta) cos(theta) 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0]
)

# Intrinsic Vector Basis that is an instance of VectorBasis
    psi = 0.0 # Define psi (angle between the x-axis and the x-axis of the Intrinsic Vector Basis)

IntrinsicVectorBasis = VectorBasis(
    Float64[cos(psi) sin(psi) 0.0;
        -sin(psi) cos(psi) 0.0;
        0.0 0.0 1.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0]
)


# Vector in 3D space defined with respect to a VectorBasis
mutable struct VectorWithBasis
    basis::VectorBasis
    v::Vector{Real}
end

# Define a constructor that takes in a VectorBasis and an optional vector argument v
function VectorWithBasisConstructor(basis::VectorBasis, v::Vector{Real}=[0.0, 0.0, 0.0])
    # Check that the input vectors have the correct dimensions
    @assert length(v) == 3 "Vector must have length 3"
    
    # Construct the VectorWithBasis object and return it
    return VectorWithBasis(basis, v)
end

# Define a default constructor that initializes the basis to the Canonical Vector Basis and the vector to zero
VectorWithBasis() = VectorWithBasis(CanonicalVectorBasis, [0.0, 0.0, 0.0])
#DefaultVectorWithBasis = VectorWithBasis()


##############################################################################################################
# SOME MODIFICATIONS TO BREAK THE CIRCULAR DEPENDENCY OF Point and ReferenceFrame

# Abstract type for reference frame
abstract type AbstractReferenceFrame end

# Point in 3D space defined with a position, velocity and acceleration and a reference to a ReferenceFrame
mutable struct Point{AbstractReferenceFrame}
    position::Vector{Real} # position vector
    velocity::Vector{Real} # velocity vector
    acceleration::Vector{Real} # acceleration vector
    reference_frame::AbstractReferenceFrame # reference frame of the point
end

# Define a constructor for Point that takes in the position, velocity, acceleration, and reference frame
function Point{AbstractReferenceFrame}(position::Vector{Real}, velocity::Vector{Real}, acceleration::Vector{Real}, reference_frame::AbstractReferenceFrame)
    return Point{AbstractReferenceFrame}(position, velocity, acceleration, reference_frame)
end

# Define a default constructor for Point that initializes the position, velocity, and acceleration to zero, and the reference frame to the inertial reference frame
function Point{AbstractReferenceFrame}()
    position = [0.0, 0.0, 0.0]
    velocity = [0.0, 0.0, 0.0]
    acceleration = [0.0, 0.0, 0.0]
    reference_frame = ReferenceFrame{AbstractReferenceFrame}()
    return Point{AbstractReferenceFrame}(position, velocity, acceleration, reference_frame)
end

# Reference Frame Type defined by a VectorBasis and an origin Point
mutable struct ReferenceFrame{AbstractReferenceFrame}
    basis::VectorBasis # basis of the reference frame
    origin::Point{AbstractReferenceFrame} # origin point of the reference frame
end

# Define a constructor for ReferenceFrame that takes in the basis and origin point
function ReferenceFrame{AbstractReferenceFrame}(basis::VectorBasis, origin::Point{AbstractReferenceFrame})
    return ReferenceFrame{AbstractReferenceFrame}(basis, origin)
end

# Define a default constructor for ReferenceFrame that initializes the basis to the canonical vector basis, and the origin to the origin point of the inertial reference frame
function ReferenceFrame{AbstractReferenceFrame}()
    basis = VectorBasis()
    origin = Point{AbstractReferenceFrame}()
    return ReferenceFrame{AbstractReferenceFrame}(basis, origin)
end

##############################################################################################################


# Particle Type defined by a Point and a mass
mutable struct Particle 
    point::Point
    mass::Real
end

# Define a constructor that takes in a Point and a mass
function ParticleConstructor(point::Point, mass::Real)
    return Particle(point, mass)
end

# Define a default constructor that initializes the point to the origin point of the inertial reference frame and the mass to zero
Particle() = Particle(Point(), 0.0) # massless particle


end # module