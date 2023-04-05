# Tests related to LinearAlgebraTypes.jl

let
    # Test MyVector object creation
    v1 = MyVector(1, 2, 3)
    @test v1.x == 1
    @test v1.y == 2
    @test v1.z == 3

    # Test MyVector object creation with invalid input
    @test_throws AssertionError MyVector(1, 2, "3")

    # Test MyVector object creation from basis
    @test MyBasis() == CanonicalBasis
    b0 = CanonicalBasis
    v2 = MyVector(1, 2, 3, b0)
    @test v2.x == 1
    @test v2.y == 2
    @test v2.z == 3

    # Test MyBasis object creation
    v_i = MyVector(1, 0, 0)
    v_j = MyVector(0, 1, 0)
    v_k = MyVector(0, 0, 1)
    omega = MyVector(0, 0, 0)
    alpha = MyVector(0, 0, 0)
    b1 = MyBasis(v_i, v_j, v_k, omega, alpha)
    @test b1.i == v_i
    @test b1.j == v_j
    @test b1.k == v_k
    @test b1.omega == omega
    @test b1.alpha == alpha

    # Test MyBasis object creation with invalid input
    @test_throws MethodError MyBasis(v_i, v_j, v_k, omega, "alpha")

    # Test MyBasis object creation from basis
    v_i = [1,0,0] # in basis b1
    v_j = [0,1,0] # in basis b1
    v_k = [0,0,1] # in basis b1
    omega = [0,0,0] # in basis b1
    alpha = [0,0,0] # in basis b1
    b2 = MyBasis(v_i, v_j, v_k, omega, alpha, b1)
    v_i_0 = MyVector(1, 0, 0, b1)
    v_j_0 = MyVector(0, 1, 0, b1)
    v_k_0 = MyVector(0, 0, 1, b1)
    omega_0 = MyVector(0, 0, 0, b1)
    alpha_0 = MyVector(0, 0, 0, b1)
    @test b2.i == v_i_0
    @test b2.j == v_j_0
    @test b2.k == v_k_0
    @test b2.omega == omega_0
    @test b2.alpha == alpha_0
end

let
    # Test MyBasis performs a proper transformation with an example
    i1 = [1/sqrt(2), 0.0, 1/sqrt(2)] # defined wrt B0
    j1 = [0.0, 1.0, 0.0] # defined wrt B0
    k1 = [-1/sqrt(2), 0.0, 1/sqrt(2)] # defined wrt B0
    omega1 = [0.0,1.0,0.0] # defined wrt B0
    alpha1 = [0.0,0.0,0.0] # defined wrt B0
    basis1 = [1/sqrt(2) 0.0 1/sqrt(2); 0.0 1.0 0.0; -1/sqrt(2) 0.0 1/sqrt(2)]
    T = inv(basis1) # transformation matrix from B1 to B0

    i2_1 = [1/sqrt(2), 0.0, 1/sqrt(2)] # defined wrt B1
    j2_1 = [0.0, 1.0, 0.0] # defined wrt B1
    k2_1 = [-1/sqrt(2), 0.0, 1/sqrt(2)] # defined wrt B1
    omega2_1 = [0.0,3.0,0.0] # defined wrt B1
    alpha2_1 = [0.0,0.0,0.0] # defined wrt B1
    #=
    i2_1 = [0.0, 0.0, 1.0] # defined wrt B1
    j2_1 = [0.0, 1.0, 0.0] # defined wrt B1
    k2_1 = [1.0, 0.0, 0.0] # defined wrt B1
    omega2_1 = [0.0,3.0,0.0] # defined wrt B1
    alpha2_1 = [0.0,0.0,0.0] # defined wrt B1
    =#

    i2_0 = T * i2_1 # defined wrt B0
    j2_0 = T * j2_1 # defined wrt B0
    k2_0 = T * k2_1 # defined wrt B0
    omega2_0 = T * omega2_1 # defined wrt B0
    alpha2_0 = T * alpha2_1 # defined wrt B0

    B1 = MyBasis(i1, j1, k1, omega1, alpha1,CanonicalBasis)
    #B1= MyBasis(MyVector(i1[1],i1[2],i1[3]), MyVector(j1[1],j1[2],j1[3]), MyVector(k1[1],k1[2],k1[3]), MyVector(omega1[1],omega1[2],omega1[3]), MyVector(alpha1[1],alpha1[2],alpha1[3]))
    #@test B1 == B1_tryingsomething
    B2 = MyBasis(i2_1, j2_1, k2_1, omega2_1, alpha2_1, B1)
    #@test B2.i.x ≈ i2_0[1] atol = 1e-8
    tol = 1e-8
    @test isapprox(B2.i.x, i2_0[1], atol=tol) 
    @test isapprox(B2.i.y, i2_0[2], atol=tol)
    @test isapprox(B2.i.z, i2_0[3], atol=tol)
    @test isapprox(B2.j.x, j2_0[1], atol=tol)
    @test isapprox(B2.j.y, j2_0[2], atol=tol)
    @test isapprox(B2.j.z, j2_0[3], atol=tol)
    @test isapprox(B2.k.x, k2_0[1], atol=tol)
    @test isapprox(B2.k.y, k2_0[2], atol=tol)
    @test isapprox(B2.k.z, k2_0[3], atol=tol)
    @test omega2_0[1] == B2.omega.x && omega2_0[2] == B2.omega.y && omega2_0[3] == B2.omega.z
    @test alpha2_0[1] == B2.alpha.x && alpha2_0[2] == B2.alpha.y && alpha2_0[3] == B2.alpha.z
    @test B2.i == MyVector(i2_1[1], i2_1[2], i2_1[3],B1)
    @test B2.j == MyVector(j2_1[1], j2_1[2], j2_1[3],B1)
    @test B2.k == MyVector(k2_1[1], k2_1[2], k2_1[3],B1)
    @test B2.omega == MyVector(omega2_1[1],omega2_1[2],omega2_1[3],B1)
    @test B2.alpha == MyVector(alpha2_1[1],alpha2_1[2],alpha2_1[3],B1)
end