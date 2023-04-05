module LinearAlgebraTypes

using LinearAlgebra

export MyVector, MyBasis, CanonicalBasis

# Define Types

# MyVector is a vector type that stores a vector in Canonical Basis. 
struct MyVector
    x
    y
    z

    # Define a constructor that takes in a Vector in canonical basis and returns a new MyVector object
    function MyVector(x0,y0,z0)
        # Input
            # x0, y0, z0 are the components of the input vector in canonical basis
        # Ouput
            # x, y, z are the components of the output vector in the canonical basis

        # Check that the input variables are numbers
        @assert typeof(x0) <: Number && typeof(y0) <: Number && typeof(z0) <: Number "Input must be a number"

        # Construct the MyVector object and return it
        x = x0
        y = y0
        z = z0
        return new(x, y, z)
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
        x = x1 * B1.i.x + y1*B1.j.x + z1*B1.k.x
        y = x1 * B1.i.y + y1 * B1.j.y + z1 * B1.k.y
        z = x1 * B1.i.z + y1 * B1.j.z + z1 * B1.k.z

        # Construct the MyVector object and return it
        return new(x, y, z)
    end
end

MyVector() = MyVector(0,0,0)


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
        basis = [i1_0.x i1_0.y i1_0.z; j1_0.x j1_0.y j1_0.z; k1_0.x k1_0.y k1_0.z]
        @assert isapprox(dot(basis[1,:],basis[2,:]), 0.0, atol=eps()) "i1_0 and j1_0 must be orthogonal"
        @assert isapprox(dot(basis[1,:],basis[3,:]), 0.0, atol=eps()) "i1_0 and k1_0 must be orthogonal"
        @assert isapprox(dot(basis[2,:],basis[3,:]), 0.0, atol=eps()) "j1_0 and k1_0 must be orthogonal"

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
        # Create MyVector objects from the input vectors in basis B1
        i2_0 = MyVector(i2_1[1],i2_1[2],i2_1[3],B1)
        j2_0 = MyVector(j2_1[1],j2_1[2],j2_1[3],B1)
        k2_0 = MyVector(k2_1[1],k2_1[2],k2_1[3],B1)
        omega2_0 = MyVector(omega2_1[1],omega2_1[2],omega2_1[3],B1)
        alpha2_0 = MyVector(alpha2_1[1],alpha2_1[2],alpha2_1[3],B1)

        # Check that the vectors are orthogonal
        basis = [i2_0.x i2_0.y i2_0.z; j2_0.x j2_0.y j2_0.z; k2_0.x k2_0.y k2_0.z]
        @assert isapprox(dot(basis[1,:],basis[2,:]), 0.0, atol=eps()) "i2_0 and j2_0 must be orthogonal"
        @assert isapprox(dot(basis[1,:],basis[3,:]), 0.0, atol=eps()) "i2_0 and k2_0 must be orthogonal"
        @assert isapprox(dot(basis[2,:],basis[3,:]), 0.0, atol=eps()) "j2_0 and k2_0 must be orthogonal"

        # Construct the MyBasis object and return it
        return new(i2_0, j2_0, k2_0, omega2_0, alpha2_0)
    end
end

i = MyVector(1.0,0.0,0.0)
println(i.x)
j = MyVector(0.0,1.0,0.0)
k = MyVector(0.0,0.0,1.0)
omega_0 = MyVector(0.0,0.0,0.0)
alpha_0 = MyVector(0.0,0.0,0.0)

MyBasis() = MyBasis(i, j, k, omega_0, alpha_0)
CanonicalBasis = MyBasis()


end # module
