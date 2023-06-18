# LinearAlgebraTypes file

using LinearAlgebra
using Plots

export MyVector, MyBasis, CanonicalBasis, i0, j0, k0, MyPoint, CanonicalOriginPoint, MyReferenceFrame, CanonicalReferenceFrame, 
        MyParticle, componentsInBasis, pos_vel_acc_inRF, plot_MyBasis, plot_MyReferenceFrame, rotation_matrix, 
        pure_rotation_MyBasis_x, pure_rotation_MyBasis_y, pure_rotation_MyBasis_z

# Define Types

# MyVector is a vector type that stores a vector in Canonical Basis. 
struct MyVector
    c0::Vector{Float64}

    # Define a constructor that takes in a Vector of length 3 in canonical basis and returns a new MyVector object
    function MyVector(c0)
        # Input
            # c0 is the input vector of length 3 in canonical basis
        # Ouput
            # c0 is the output vector of length 3 in canonical basis

        # Check that the input vectors are length 3
        @assert length(c0) == 3 "Input vector must be length 3"

        # Construct the MyVector object and return it
        return new(c0)
    end
    # Define a constructor that takes in a Vector of length 3 in any basis and that basis in canonical basis and returns a new MyVector object
    function MyVector(c1,B1)
        # Input
            # c1 is the input vector of length 3 in basis B1
            # B1 is a MyBasis object that represents the basis B1 in the canonical basis B0 (could also be called B1_0)
        # Ouput
            # c0 is the output vector of length 3 in the canonical basis

        # Check that the input vectors are length 3
        @assert length(c1) == 3 "Input vector must be length 3"

        # Convert the input vector from basis B1 to the canonical basis
        c0 = c1[1] * B1.i.c0 + c1[2] * B1.j.c0 + c1[3] * B1.k.c0

        # Construct the MyVector object and return it
        return new(c0)
    end

    # Define a constructor that takes in a Vector in canonical basis and returns a new MyVector object
    function MyVector(x0,y0,z0)
        # Input
            # x0, y0, z0 are the components of the input vector in canonical basis
        # Ouput
            # x, y, z are the components of the output vector in the canonical basis

        # Check that the input variables are numbers
        @assert typeof(x0) <: Number && typeof(y0) <: Number && typeof(z0) <: Number "Input must be a number"

        # Construct the MyVector object and return it
        c0 = [x0,y0,z0] 
        return new(c0)
    end
    # Define a constructor that takes in a Vector in any basis and that basis in canonical basis and returns a new MyVector object
    function MyVector(x1,y1,z1,B1)
        # Input
            # x1, y1, z1 are the components of the input vector in basis B1
            # B1 is a MyBasis object that represents the basis B1 in the canonical basis B0 (could also be called B1_0)
        # Ouput
            # x, y, z are the components of the output vector in the canonical basis
        
        # Check that the input variables are numbers
        @assert typeof(x1) <: Number && typeof(y1) <: Number && typeof(z1) <: Number "Input must be a number"

        # Convert the input vector from basis B1 to the canonical basis
        x0 = x1 * B1.i.c0[1] + y1 * B1.j.c0[1] + z1 * B1.k.c0[1]
        y0 = x1 * B1.i.c0[2] + y1 * B1.j.c0[2] + z1 * B1.k.c0[2]
        z0 = x1 * B1.i.c0[3] + y1 * B1.j.c0[3] + z1 * B1.k.c0[3]

        # Construct the MyVector object and return it
        c0 = [x0,y0,z0] 
        return new(c0)
    end
end

MyVector() = MyVector(0,0,0)
Base.:+(a::MyVector, b::MyVector) = MyVector(a.c0 + b.c0)
Base.:-(a::MyVector, b::MyVector) = MyVector(a.c0 - b.c0)
Base.:-(a::MyVector) = MyVector(-a.c0)
Base.:(==)(a::MyVector, b::MyVector) = a.c0 == b.c0
Base.:*(a::Number, b::MyVector) = MyVector(a * b.c0)
Base.:*(a::MyVector, b::Number) = MyVector(a.c0 * b)
Base.:/(a::MyVector, b::Number) = MyVector(a.c0 / b)
Base.:/(a::MyVector, b::MyVector) = MyVector(a.c0 ./ b.c0)
LinearAlgebra.:(dot)(a::MyVector, b::MyVector) = dot(a.c0,b.c0)
LinearAlgebra.:(cross)(a::MyVector, b::MyVector) = MyVector(cross(a.c0,b.c0))
LinearAlgebra.:(norm)(a::MyVector) = norm(a.c0)
#Base.:(isapprox)(a::MyVector, b::MyVector) = isapprox(a.c0,b.c0)
Base.:(isapprox)(a::MyVector, b::MyVector; atol::Real, rtol::Real) = isapprox(a.c0,b.c0;atol,rtol)
Base.:(getindex)(a::MyVector, i::Int) = a.c0[i]
Base.:(length)(a::MyVector) = length(a.c0)

# MyBasis is a basis type that stores a basis in Canonical Basis, this is its orthonormal vectors, and its angular velocity and acceleration.
struct MyBasis
    i::MyVector
    j::MyVector
    k::MyVector
    omega::MyVector
    alpha::MyVector

    # Define a constructor that takes in 5 MyVector objects in canonical basis and returns a new MyBasis object
    function MyBasis(i1_0::MyVector, j1_0::MyVector, k1_0::MyVector, omega1_0::MyVector, alpha1_0::MyVector)
        # Input
            # i, j, k are the orthonormal vectors of the basis in canonical basis
            # omega is the angular velocity of the basis in canonical basis
            # alpha is the angular acceleration of the basis in canonical basis
        # Ouput
            # i, j, k are the orthonormal vectors of the basis in canonical basis
            # omega is the angular velocity of the basis in canonical basis
            # alpha is the angular acceleration of the basis in canonical basis

        # Check that the input vectors are of the correct type
        @assert typeof(i1_0) == MyVector "i1_0 must be a MyVector object"
        @assert typeof(j1_0) == MyVector "j1_0 must be a MyVector object"
        @assert typeof(k1_0) == MyVector "k1_0 must be a MyVector object"
        @assert typeof(omega1_0) == MyVector "omega1_0 must be a MyVector object"
        @assert typeof(alpha1_0) == MyVector "alpha1_0 must be a MyVector object"
        
        # Check that the vectors are orthogonal
        @assert isapprox(dot(i1_0,j1_0), 0.0, atol=eps()) "i1_0 and j1_0 must be orthogonal"
        @assert isapprox(dot(i1_0,k1_0), 0.0, atol=eps()) "i1_0 and k1_0 must be orthogonal"
        @assert isapprox(dot(j1_0,k1_0), 0.0, atol=eps()) "j1_0 and k1_0 must be orthogonal"

        # Construct the MyBasis object and return it
        return new(i1_0, j1_0, k1_0, omega1_0, alpha1_0)
    end

    # Define a constructor that takes in 5 Vectors in any basis and that basis in canonical basis and returns a new MyBasis object
    function MyBasis(i2_1::Vector, j2_1::Vector, k2_1::Vector, omega2_1::Vector, alpha2_1::Vector, B1::MyBasis)
        #Input
            # i2_1, j2_1, k2_1 are the orthonormal vectors of the basis B2 in basis B1
            # omega2_1 is the angular velocity of the basis B2 in basis B1
            # alpha2_1 is the angular acceleration of the basis B2 in basis B1
            # B1 is a MyBasis object that represents the basis B1 in the canonical basis B0
        # Ouput
            # i, j, k are the orthonormal vectors of the basis B2 in canonical basis B0
            # omega is the angular velocity of the basis B2 in canonical basis B0
            # alpha is the angular acceleration of the basis B2 in canonical basis B0

        # Check that the input vectors are of the correct length and type
        @assert typeof(B1) == MyBasis "B1 must be a MyBasis object"
        @assert length(i2_1) == 3 "i2_1 must have length 3"
        @assert length(j2_1) == 3 "j2_1 must have length 3"
        @assert length(k2_1) == 3 "k2_1 must have length 3"
        @assert length(omega2_1) == 3 "omega2_1 must have length 3"
        @assert length(alpha2_1) == 3 "alpha2_1 must have length 3"

        # Convert the input basis B2 from basis B1 to the canonical basis
        # Create MyVector objects from the input basis vectors in basis B1
        i2_0 = MyVector(i2_1[1],i2_1[2],i2_1[3],B1)
        j2_0 = MyVector(j2_1[1],j2_1[2],j2_1[3],B1)
        k2_0 = MyVector(k2_1[1],k2_1[2],k2_1[3],B1)

        # Convert omega2_1 to MyVector before sum
        omega2_1 = MyVector(omega2_1,B1)

        # Convert alpha2_1 to MyVector before sum
        alpha2_1 = MyVector(alpha2_1,B1)

        # Convert the input angular velocity and acceleration from basis B1 to the canonical basis
        omega2_0 = omega2_1 + B1.omega
        alpha2_0 = alpha2_1 + cross(B1.omega, omega2_1) + B1.alpha

        # Check that the vectors are orthogonal
        @assert isapprox(dot(i2_0,j2_0), 0.0, atol=eps()) "i2_0 and j2_0 must be orthogonal"
        @assert isapprox(dot(i2_0,k2_0), 0.0, atol=eps()) "i2_0 and k2_0 must be orthogonal"
        @assert isapprox(dot(j2_0,k2_0), 0.0, atol=eps()) "j2_0 and k2_0 must be orthogonal"

        # Construct the MyBasis object and return it
        return new(i2_0, j2_0, k2_0, omega2_0, alpha2_0)
    end
end

i0 = MyVector(1.0,0.0,0.0)
j0 = MyVector(0.0,1.0,0.0)
k0 = MyVector(0.0,0.0,1.0)
omega_0 = MyVector(0.0,0.0,0.0)
alpha_0 = MyVector(0.0,0.0,0.0)

MyBasis() = MyBasis(i0, j0, k0, omega_0, alpha_0)
CanonicalBasis = MyBasis()


struct MyPoint 
    position::MyVector
    velocity::MyVector
    acceleration::MyVector

    # Define a constructor that takes in 3 MyVector objects (in canonical basis) and returns a MyPoint object
    function MyPoint(pos::MyVector,veloc::MyVector,accel::MyVector)
        #Input
            # pos is the position of the point in canonical basis B0
            # veloc is the velocity of the point in canonical basis B0
            # accel is the acceleration of the point in canonical basis B0
        # Ouput
            # position is the position of the point in canonical basis B0
            # velocity is the velocity of the point in canonical basis B0
            # acceleration is the acceleration of the point in canonical basis B0

        # Check that the input vectors are of the correct length and type
        @assert typeof(pos) == MyVector "pos must be a MyVector object"
        @assert typeof(veloc) == MyVector "veloc must be a MyVector object"
        @assert typeof(accel) == MyVector "accel must be a MyVector object"

        # Construct the MyPoint object and return it
        return new(pos, veloc, accel)
    end   
end

MyPoint() = MyPoint(MyVector(0.0,0.0,0.0),MyVector(0.0,0.0,0.0),MyVector(0.0,0.0,0.0))

CanonicalOriginPoint = MyPoint()


struct MyReferenceFrame 
    basis::MyBasis
    origin::MyPoint

    # Define a constructor that takes in a MyBasis object and a MyPoint object and returns a MyReferenceFrame object
    function MyReferenceFrame(basis::MyBasis,origin::MyPoint)
        #Input
            # basis is a MyBasis object that represents the basis of the reference frame in the canonical basis B0
            # origin is a MyPoint object that represents the origin of the reference frame in the canonical basis B0
        # Ouput
            # basis is a MyBasis object that represents the basis of the reference frame in the canonical basis B0
            # origin is a MyPoint object that represents the origin of the reference frame in the canonical basis B0

        # Check that the input vectors are of the correct length and type
        @assert typeof(basis) == MyBasis "basis must be a MyBasis object"
        @assert typeof(origin) == MyPoint "origin must be a MyPoint object"

        # Construct the MyReferenceFrame object and return it
        return new(basis, origin)
    end 
end

MyReferenceFrame() = MyReferenceFrame(CanonicalBasis,CanonicalOriginPoint)
CanonicalReferenceFrame = MyReferenceFrame()


# Define a constructor that takes in 3 Vector objects and the Reference Frame they are in and returns a MyPoint object
function MyPoint(pos1::Vector,veloc1::Vector,accel1::Vector,RF1::MyReferenceFrame)
    #Input
        # pos1 is the position of the point in reference frame R1
        # veloc1 is the velocity of the point in reference frame R1
        # accel1 is the acceleration of the point in reference frame R1
        # R1 is a MyReferenceFrame object that represents the reference frame R1 in the canonical basis B0
    # Ouput
        # position is the position of the point in canonical basis B0
        # velocity is the velocity of the point in canonical basis B0
        # acceleration is the acceleration of the point in canonical basis B0

    # Check that the input vectors are of the correct length and type
    @assert typeof(RF1) == MyReferenceFrame "R1 must be a MyReferenceFrame object"
    @assert length(pos1) == 3 "pos1 must have length 3"
    @assert length(veloc1) == 3 "veloc1 must have length 3"
    @assert length(accel1) == 3 "accel must have length 3"

    # Convert the input vectors to MyVector objects so that they are all represented in canonical components
    pos1 = MyVector(pos1,RF1.basis)
    veloc1 = MyVector(veloc1,RF1.basis)
    accel1 = MyVector(accel1,RF1.basis)

    # Convert the input vectors from reference frame RF1 to the Canonical Reference Frame
    pos0 = RF1.origin.position + pos1
        # Cross product of omega10 x pos1
        cross1 = cross(RF1.basis.omega,pos1)
    veloc0 = RF1.origin.velocity + veloc1 + cross1
        # Cross product of alpha10 x pos10
        cross2 = cross(RF1.basis.alpha,pos1)
        # Cross product of omega10 x cross1
        cross3 = cross(RF1.basis.omega,cross1)
        # Cross product of omega10 x veloc1
        cross4 = cross(RF1.basis.omega,veloc1)
    accel0 = RF1.origin.acceleration + accel1 + cross2 + cross3 + 2*cross4

    # Construct the MyPoint object and return it
    return MyPoint(pos0, veloc0, accel0)
end


struct MyParticle
    point::MyPoint
    mass::Number

    # Define a constructor that takes in a MyPoint objects and the mass and returns a MyParticle object
    function MyParticle(point::MyPoint,mass::Number)
        #Input 
            # point is a MyPoint object that represents the position, velocity, and acceleration of the particle
            # mass is the mass of the particle
        # Ouput
            # point is a MyPoint object that represents the position, velocity, and acceleration of the particle
            # mass is the mass of the particle
        
        # Check that the inputs are of the correct type
        @assert typeof(point) == MyPoint "point must be a MyPoint object"
        @assert typeof(mass) <: Number "mass must be a number"

        # Construct the MyParticle object and return it
        return new(point, mass)
    end

end


# Function that given a MyVector object projects it into another base
function componentsInBasis(v::MyVector,B1::MyBasis)
    #Input
        # v is a MyVector object that represents the vector in the canonical basis B0
        # B1 is a MyBasis object that represents the basis in which the vector is to be projected
    # Ouput
        # v1 is the vector that represents the vector v in the basis B1

    # Check that the input vectors are of the correct length and type
    @assert typeof(v) == MyVector "v must be a MyVector object"
    @assert typeof(B1) == MyBasis "B must be a MyBasis object"

    # Project the vector into the new basis
    # v1 = [dot(v.c0,B1.i.c0),dot(v.c0,B1.j.c0),dot(v.c0,B1.k.c0)] # this would only work for orthonormal 
    transform_matrix = inv([B1.i[1] B1.j[1] B1.k[1]; B1.i[2] B1.j[2] B1.k[2]; B1.i[3] B1.j[3] B1.k[3]])
    # Project the vector into the new basis
    v1 = transform_matrix*[v[1] v[2] v[3]]'
    v1 = [v1[1],v1[2],v1[3]]

    # Return the projected vector
    return v1
end

# Function that given a MyPoint object gives the Point position, velocity and acceleration in a given Reference Frame (RF1)
function pos_vel_acc_inRF(point::MyPoint,RF1::MyReferenceFrame)
    #Input
        # point is a MyPoint object that represents the position, velocity, and acceleration of the point
        # RF1 is a MyReferenceFrame object that represents the reference frame in which the point is to be projected
    # Ouput
        # pos1 is the position of the point in reference frame RF1
        # veloc1 is the velocity of the point in reference frame RF1
        # accel1 is the acceleration of the point in reference frame RF1

    # Check that the input vectors are of the correct length and type
    @assert typeof(point) == MyPoint "point must be a MyPoint object"
    @assert typeof(RF1) == MyReferenceFrame "RF1 must be a MyReferenceFrame object"

    # Convert the point from the canonical reference frame to the reference frame RF1
    pos10 = point.position - RF1.origin.position
        # Cross product of omega10 x pos1
        cross1 = cross(RF1.basis.omega,pos10)
    veloc10 = point.velocity - RF1.origin.velocity - cross1
        # Cross product of alpha10 x pos10
        cross2 = cross(RF1.basis.alpha,pos10)
        # Cross product of omega10 x cross1
        cross3 = cross(RF1.basis.omega,cross1)
        # Cross product of omega10 x veloc1
        cross4 = cross(RF1.basis.omega,veloc10)
    accel10 = point.acceleration - RF1.origin.acceleration - cross2 - cross3 - 2*cross4

    pos1 = componentsInBasis(pos10,RF1.basis)
    veloc1 = componentsInBasis(veloc10,RF1.basis)
    accel1 = componentsInBasis(accel10,RF1.basis)

    # Return the point in the reference frame RF1
    return pos1, veloc1, accel1
end

# Function that given a MyBasis, plots it in 3D using Plots.jl
function plot_MyBasis(basis::MyBasis)
    i = basis.i.c0
    j = basis.j.c0
    k = basis.k.c0

    p = plot3d(xlim=(-1, 1), ylim=(-1, 1), zlim=(-1, 1), aspect_ratio=:equal)
    plot!(p, [0, i[1]], [0, i[2]], [0, i[3]], arrow=true, label="i", linecolor=:red, arrowcolor=:red)
    plot!(p, [0, j[1]], [0, j[2]], [0, j[3]], arrow=true, label="j", linecolor=:green, arrowcolor=:green)
    plot!(p, [0, k[1]], [0, k[2]], [0, k[3]], arrow=true, label="k", linecolor=:blue, arrowcolor=:blue)

    xlabel!("X")
    ylabel!("Y")
    zlabel!("Z")
    title!("Basis Vectors")
    
    display(p)
end

# Recipe to plot a MyBasis just by using plot(b::MyBasis)
@recipe function plot(b::MyBasis)
    i = b.i.c0
    j = b.j.c0
    k = b.k.c0

    p = plot3d!(xlim=(-1, 1), ylim=(-1, 1), zlim=(-1, 1), aspect_ratio=:equal)
    plot!(p, [0, i[1]], [0, i[2]], [0, i[3]], arrow=true, label="i", linecolor=:red, arrowcolor=:red)
    plot!(p, [0, j[1]], [0, j[2]], [0, j[3]], arrow=true, label="j", linecolor=:green, arrowcolor=:green)
    plot!(p, [0, k[1]], [0, k[2]], [0, k[3]], arrow=true, label="k", linecolor=:blue, arrowcolor=:blue)

    xlabel!("X")
    ylabel!("Y")
    zlabel!("Z")
    title!("Basis Vectors")
    
    display(p)
end

# Recipe to plot a MyReferenceFrame just by using plot(RF::MyReferenceFrame)
@recipe function plot(f::MyReferenceFrame)
    i = f.basis.i.c0
    j = f.basis.j.c0
    k = f.basis.k.c0

    origin = f.origin.position

    p = plot3d(aspect_ratio=:equal)


    plot!(p, [origin[1], origin[1] + i[1]], [origin[2], origin[2] + i[2]], [origin[3], origin[3] + i[3]], arrow=true, label="i", linecolor=:red, arrowcolor=:red)
    plot!(p, [origin[1], origin[1] + j[1]], [origin[2], origin[2] + j[2]], [origin[3], origin[3] + j[3]], arrow=true, label="j", linecolor=:green, arrowcolor=:green)
    plot!(p, [origin[1], origin[1] + k[1]], [origin[2], origin[2] + k[2]], [origin[3], origin[3] + k[3]], arrow=true, label="k", linecolor=:blue, arrowcolor=:blue)
    plot!(p, [origin[1]], [origin[2]], [origin[3]], marker=:circle, markersize=5, label="Origin")
    
    xlabel!("X")
    ylabel!("Y")
    zlabel!("Z")
    title!("Reference Frame")
    
    display(p)
end

# Function to export the recipe outside of the module
function plot_MyReferenceFrame(f::MyReferenceFrame)
    plot(f)
end

# Function that given a MyBasis, gives the rotation transform_matrix
function rotation_matrix(basis::MyBasis)
    i = basis.i.c0
    j = basis.j.c0
    k = basis.k.c0

    return R = [i[1] j[1] k[1]; i[2] j[2] k[2]; i[3] j[3] k[3]]
end

# Function that creates a new MyBasis object which is a pure rotation of the input MyBasis object with respect to the x-axis given the angle of rotation
function pure_rotation_MyBasis_x(basis::MyBasis,theta::Number)
    # Input:
        # basis is a MyBasis object that represents the basis to be rotated
        # theta is the angle of rotation in radians
    # Output:
        # basis is a MyBasis object that represents the basis to be rotated

    # Get the rotation transform_matrix
    R = rotation_matrix(basis)

    # Rotate the basis with respect to the x-axis
        i = R*[1.0,0.0,0.0]
        j = R*[0.0,cos(theta),sin(theta)]
        k = R*[0.0,-sin(theta),cos(theta)]

    # Construct the MyBasis object and return it
    return MyBasis(MyVector(i),MyVector(j),MyVector(k),basis.omega,basis.alpha)
end

# Function that creates a new MyBasis object which is a pure rotation of the input MyBasis object with respect to the y-axis given the angle of rotation
function pure_rotation_MyBasis_y(basis::MyBasis,theta::Number)
    # Input:
        # basis is a MyBasis object that represents the basis to be rotated
        # theta is the angle of rotation in radians
    # Output:
        # basis is a MyBasis object that represents the basis to be rotated

    # Get the rotation transform_matrix
    R = rotation_matrix(basis)

    # Rotate the basis with respect to the y-axis
        i = R*[cos(theta),0.0,-sin(theta)]
        j = R*[0.0,1.0,0.0]
        k = R*[sin(theta),0.0,cos(theta)]

    # Construct the MyBasis object and return it
    return MyBasis(MyVector(i),MyVector(j),MyVector(k),basis.omega,basis.alpha)
end

# Function that creates a new MyBasis object which is a pure rotation of the input MyBasis object with respect to the z-axis given the angle of rotation
function pure_rotation_MyBasis_z(basis::MyBasis,theta::Number)
    # Input:
        # basis is a MyBasis object that represents the basis to be rotated
        # theta is the angle of rotation in radians
    # Output:
        # basis is a MyBasis object that represents the basis to be rotated

    # Get the rotation transform_matrix
    R = rotation_matrix(basis)

    # Rotate the basis with respect to the z-axis
        i = R*[cos(theta),sin(theta),0.0]
        j = R*[-sin(theta),cos(theta),0.0]
        k = R*[0.0,0.0,1.0]

    # Construct the MyBasis object and return it
    return MyBasis(MyVector(i),MyVector(j),MyVector(k),basis.omega,basis.alpha)
end

