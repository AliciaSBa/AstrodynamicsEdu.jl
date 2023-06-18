# LINEAR ALGEBRA TYPES TEST

using AstrodynamicsEdu.LinearAlgebraTypes
using LinearAlgebra
using Test

# Test MyVector constructors
function test_MyVector()
    # Test the constructor that takes in a vector in canonical basis
    v1 = MyVector(1, 2, 3)
    @test v1.c0 == [1, 2, 3]
    v2 = MyVector([1, 2, 3])
    @test v2.c0 == [1, 2, 3]

    # Test MyVector object creation with invalid input
    @test_throws AssertionError MyVector(1, 2, "3")
    @test_throws AssertionError MyVector([1,2,3,4])

    # Test the constructor that takes in a vector in any basis
    B1 = MyBasis(MyVector([1, 0, 0]), MyVector([0, 1, 0]), MyVector([0, 0, 1]), MyVector([0, 0, 0]), MyVector([0, 0, 0]))
    v2 = MyVector([1, 2, 3], B1)
    @test v2.c0 == [1, 2, 3]

    # Test the constructor that takes in a vector in any basis and converts it to canonical basis
    v3 = MyVector(1, 2, 3, B1)
    @test v3.c0 == [1, 2, 3]
    v4 = MyVector([1, 2, 3], B1)
    @test v4.c0 == [1, 2, 3]

    # Test MyVector object creation from basis
    @test MyBasis() == CanonicalBasis
    b0 = CanonicalBasis
    v5 = MyVector(1, 2, 3, b0)
    @test v5.c0 == [1, 2, 3]

    # Test the MyVector default constructor
    v = MyVector()
    @test v.c0 == [0, 0, 0]

end

function test_MyVector_operations()
    @test MyVector(1.0, 2.0, 3.0) + MyVector(1.0, 2.0, 3.0) == MyVector(2.0, 4.0, 6.0)
    @test MyVector(1.0, 2.0, 3.0) - MyVector(1.0, 2.0, 3.0) == MyVector(0.0, 0.0, 0.0)
    @test -MyVector(1.0, 2.0, 3.0) == MyVector(-1.0, -2.0, -3.0)
    @test MyVector(1.0, 2.0, 3.0) * 2.0 == MyVector(2.0, 4.0, 6.0)
    @test 2.0 * MyVector(1.0, 2.0, 3.0) == MyVector(2.0, 4.0, 6.0)
    @test MyVector(1.0, 2.0, 3.0) * 2.0 == MyVector(2.0, 4.0, 6.0)
    @test MyVector(1.0, 2.0, 3.0) / 2.0 == MyVector(0.5, 1.0, 1.5)
    @test dot(MyVector(1.0, 2.0, 3.0), MyVector(1.0, 2.0, 3.0)) == 14.0
    @test cross(MyVector(1.0, 2.0, 3.0), MyVector(1.0, 2.0, 3.0)) == MyVector(0.0, 0.0, 0.0)
    @test norm(MyVector(1.0, 2.0, 3.0)) == sqrt(14.0)
    #@test isapprox(MyVector(1.0, 2.0, 3.0), MyVector(1.0, 2.0, 3.0)) == true
    #@test isapprox(MyVector(1.0, 2.0, 3.0), MyVector(1.0, 2.0, 3.1)) == false
    @test isapprox(MyVector(1.0, 2.0, 3.0), MyVector(1.0, 2.0, 3.1); atol=0.2, rtol=0.05) == true
    @test getindex(MyVector(1.0, 2.0, 3.0), 1) == 1.0
    @test length(MyVector(1.0, 2.0, 3.0)) == 3
end


# Test MyBasis constructors
function test_MyBasis()
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

    # Test construction of a MyBasis object in the canonical basis
    b0 = CanonicalBasis
    i0 = MyVector(1,0,0)
    j0 = MyVector(0,1,0)
    k0 = MyVector(0,0,1)
    omega0 = MyVector(0,0,0)
    alpha0 = MyVector(0,0,0)
    b00 = MyBasis(i0,j0,k0,omega0,alpha0)
    @test b00.i == b0.i
    @test b00.j == b0.j
    @test b00.k == b0.k
    @test b00.omega == omega0
    @test b00.alpha == alpha0
    
    # Test construction of a MyBasis object in a non-canonical basis
    i_21 = [1,0,0] # in basis b1
    j_21 = [0,1,0] # in basis b1
    k_21 = [0,0,1] # in basis b1
    omega_21 = [0, 2, 0] # in basis b1
    alpha_21 = [0, 0, 1] # in basis b1
    b2 = MyBasis(i_21, j_21, k_21, omega_21, alpha_21, b1)
    i_20 = MyVector(1, 0, 0, b1)
    j_20 = MyVector(0, 1, 0, b1)
    k_20 = MyVector(0, 0, 1, b1)
    omega_21_0 = MyVector(omega_21,b1)
    omega_20 = omega_21_0 + b1.omega
    alpha_21_0 = MyVector(alpha_21,b1)
    alpha_20 = alpha_21_0 + cross(b1.omega, omega_21_0) + b1.alpha
    @test b2.i == i_20
    @test b2.j == j_20
    @test b2.k == k_20
    @test b2.omega == omega_20
    @test b2.alpha == alpha_20

    # Test MyBasis performs a proper transformation with an Example
 
end

# Test MyPoint constructors
function test_MyPoint()
    # Test MyPoint object creation
    pos10 = MyVector(1,1,1)
    vel10 = MyVector(10,0,0)
    acc10 = MyVector(-2,0,0)
    point1 = MyPoint(pos10, vel10, acc10)
    @test point1.position == pos10
    @test point1.velocity == vel10
    @test point1.acceleration == acc10

    # Test MyPoint object creation with invalid input
    @test_throws MethodError MyPoint(1,2,3)
    @test_throws MethodError MyPoint(pos10, vel10, "acc10")
    @test_throws MethodError MyPoint(pos10, vel10, [-2,0,0])

    # Test the constructor that takes in a point with a reference frame
    pos20 = [1,1,1]
    vel20 = [10,0,0]
    acc20 = [-2,0,0]
    point2 = MyPoint(pos20,vel20,acc20,CanonicalReferenceFrame)
    @test point2.position.c0 == pos20
    @test point2.velocity.c0 == vel20
    @test point2.acceleration.c0 == acc20

    # Test the CanonicalOriginPoint
    emptyPoint = MyPoint()
    @test emptyPoint.position == CanonicalOriginPoint.position
    @test emptyPoint.velocity == CanonicalOriginPoint.velocity
    @test emptyPoint.acceleration == CanonicalOriginPoint.acceleration
end

# Test MyReferenceFrame constructors
function test_MyReferenceFrame()
    # Test MyReferenceFrame object creation
    v_i = MyVector(1, 0, 0)
    v_j = MyVector(0, 1, 0)
    v_k = MyVector(0, 0, 1)
    omega = MyVector(0, 0, 0)
    alpha = MyVector(0, 0, 0)
    b1 = MyBasis(v_i, v_j, v_k, omega, alpha)
    origin = MyPoint()
    rf1 = MyReferenceFrame(b1, origin)
    @test rf1.basis == b1
    @test rf1.origin == origin

    # Test MyReferenceFrame object creation with invalid input
    @test_throws MethodError MyReferenceFrame(b1, "origin")
    @test_throws MethodError MyReferenceFrame(b1, [0,0,0])
    @test_throws MethodError MyReferenceFrame("b1", origin)
    @test_throws MethodError MyReferenceFrame([1,0,0], origin)

    # Test the CanonicalReferenceFrame
    emptyRF = MyReferenceFrame()
    @test emptyRF.basis == CanonicalReferenceFrame.basis
    @test emptyRF.origin == CanonicalReferenceFrame.origin
end

# Test MyParticle constructors
function test_MyParticle()
    # Test MyParticle object creation
    pos10 = MyVector(1,1,1)
    vel10 = MyVector(10,0,0)
    acc10 = MyVector(-2,0,0)
    point1 = MyPoint(pos10, vel10, acc10)
    mass1 = 1 # kg
    particle1 = MyParticle(point1, mass1)
    @test particle1.point.position == pos10
    @test particle1.point.velocity == vel10
    @test particle1.point.acceleration == acc10
    @test particle1.mass == mass1

    # Test MyParticle object creation with invalid input
    @test_throws MethodError MyParticle(1,2)
    @test_throws AssertionError MyParticle(point1, "mass1")
    @test_throws AssertionError MyParticle(point1, [1])

end

# Test function componentsInBasis
function test_componentsInBasis()
    # Test componentInBasis with a MyVector
    i10 = MyVector(0, 0, 1)
    j10 = MyVector(0, 1, 0)
    k10 = MyVector(1, 0, 0)
    omega10 = MyVector(0, 2, 0)
    alpha10 = MyVector(0, 0, 1)
    b1 = MyBasis(i10, j10, k10, omega10, alpha10)
    v11 = [1, 1, 1]
    v10 = MyVector(v11, b1)
    @test componentsInBasis(v10, b1) == v11
    i21 = [0, 1, 0]
    j21 = [1, 0, 0]
    k21 = [0, 0, 1]
    omega21 = [0, 0, 2]
    alpha21 = [1, 0, 0]
    b2 = MyBasis(i21, j21, k21, omega21, alpha21,b1)
    v21 = [1, 1, 1]
    v20 = MyVector(v21, b2)
    @test componentsInBasis(v20, b1) == v21

    # Test componentInBasis with invalid input
    @test_throws MethodError componentsInBasis(v11, b1)
    @test_throws MethodError componentsInBasis(v10, "b1")
    @test_throws MethodError componentsInBasis(v10, [1,0,0])
end

# Test function pos_vel_acc_inRF
function test_pos_vel_acc_inRF()
    # Test pos_vel_acc_inRF with a MyPoint
    i10 = MyVector(0, 0, 1)
    j10 = MyVector(0, 1, 0)
    k10 = MyVector(1, 0, 0)
    omega10 = MyVector(0, 2, 0)
    alpha10 = MyVector(0, 0, 1)
    b1 = MyBasis(i10, j10, k10, omega10, alpha10)
    origin1 = MyPoint()
    rf1 = MyReferenceFrame(b1, origin1)
    pos11 = [1,1,1]
    vel11 = [10,0,0]
    acc11 = [-2,0,0]
    point10 = MyPoint(pos11, vel11, acc11,rf1)
    pos1,vel1,acc1 = pos_vel_acc_inRF(point10, rf1)
    @test pos1 == pos11
    @test vel1 == vel11
    @test acc1 == acc11

    # Test pos_vel_acc_inRF with invalid input
    @test_throws MethodError pos_vel_acc_inRF(1, rf1)
    @test_throws UndefVarError pos_vel_acc_inRF(point1, "rf1")
end

# Test function rotation_matrix
function test_rotation_matrix()
    # Test rotation_matrix with a MyBasis
    i10 = MyVector(0, 0, 1)
    j10 = MyVector(0, 1, 0)
    k10 = MyVector(1, 0, 0)
    omega10 = MyVector(0, 2, 0)
    alpha10 = MyVector(0, 0, 1)
    b1 = MyBasis(i10, j10, k10, omega10, alpha10)
    R1 = rotation_matrix(b1)
    @test R1 == [0 0 1; 0 1 0; 1 0 0]

    # Test rotation_matrix with invalid input
    @test_throws MethodError rotation_matrix(1)
    @test_throws MethodError rotation_matrix("b1")
end

# Test functions of pure rotation of MyBasis
function test_pure_rotation_MyBasis()
    i10 = MyVector(0, 0, 1)
    j10 = MyVector(0, 1, 0)
    k10 = MyVector(1, 0, 0)
    omega10 = MyVector(0, 2, 0)
    alpha10 = MyVector(0, 0, 1)
    b1 = MyBasis(i10, j10, k10, omega10, alpha10)
    # Test pure_rotation_MyBasis_x
    b1_rotated_x = pure_rotation_MyBasis_x(b1, pi/2)
    @test b1_rotated_x.i == i10
    @test b1_rotated_x.j == j10*cos(pi/2) + k10*sin(pi/2)
    @test b1_rotated_x.k == -j10*sin(pi/2) + k10*cos(pi/2)
    @test b1_rotated_x.omega == omega10
    @test b1_rotated_x.alpha == alpha10
    # Test pure_rotation_MyBasis_y
    b1_rotated_y = pure_rotation_MyBasis_y(b1, pi/2)
    @test b1_rotated_y.i == i10*cos(pi/2) - k10*sin(pi/2)
    @test b1_rotated_y.j == j10
    @test b1_rotated_y.k == i10*sin(pi/2) + k10*cos(pi/2)
    @test b1_rotated_y.omega == omega10
    @test b1_rotated_y.alpha == alpha10
    # Test pure_rotation_MyBasis_z
    b1_rotated_z = pure_rotation_MyBasis_z(b1, pi/2)
    @test b1_rotated_z.i == i10*cos(pi/2) + j10*sin(pi/2)
    @test b1_rotated_z.j == -i10*sin(pi/2) + j10*cos(pi/2)
    @test b1_rotated_z.k == k10
    @test b1_rotated_z.omega == omega10
    @test b1_rotated_z.alpha == alpha10
end


# Test LinearAlgebraTypes with an Example
function test_LinearAlgebraTypes()
    # rf0 = CanonicalReferenceFrame
    # origin0 = CanonicalOriginPoint

    # b1 = Basis of rf1, known wrt b0
    i10 = MyVector(0, 0, 1)
    j10 = MyVector(0, 1, 0)
    k10 = MyVector(1, 0, 0)
    omega10 = MyVector(0, 2, 0)
    alpha10 = MyVector(0, 0, 1)
    b1 = MyBasis(i10, j10, k10, omega10, alpha10)

    # Origin Point of rf1 wrt rf0
    posO1_10 = MyVector(0,3,1)
    velO1_10 = MyVector(0,0,5)
    accO1_10 = MyVector(0,0,0)
    origin1 = MyPoint(posO1_10,velO1_10,accO1_10) # wrt Canonical Reference Frame
    
    # rf1 = ReferenceFrame wrt rf0
    rf1 = MyReferenceFrame(b1, origin1)

    # p1 = Point wrt rf1
    posP1_11 = [1,1,1]
    velP1_11 = [10,0,0]
    accP1_11 = [-2,0,0]
    p1 = MyPoint(posP1_11, velP1_11, accP1_11, rf1)

    # particle1 = Particle at point 1 wrt rf1
    mass1 = 10 # kg
    particle1 = MyParticle(p1, mass1)

    # b2 = Basis of rf2, known wrt b1
    i21 = [0, 1, 0]
    j21 = [1, 0, 0]
    k21 = [0, 0, 1]
    omega21 = [0, 0, 3]
    alpha21 = [0, 0, 1]
    b2 = MyBasis(i21, j21, k21, omega21, alpha21, b1)

    # Origin Point of rf2 wrt rf1
    posO2_21 = [-4,0,0]
    velO2_21 = [0,0,10]
    accO2_21 = [0,0,1]
    origin2 = MyPoint(posO2_21,velO2_21,accO2_21,rf1) # wrt rf1

    # rf2 = ReferenceFrame wrt rf1
    rf2 = MyReferenceFrame(b2, origin2)

    # p2 = Point wrt rf2
    posP2_22 = [1,1,1]
    velP2_22 = [20,0,0]
    accP2_22 = [-2,0,0]
    p2 = MyPoint(posP2_22, velP2_22, accP2_22, rf2)

    # particle2 = Particle at point 2 wrt rf2
    mass2 = 5 # kg
    particle2 = MyParticle(p2, mass2)

    # Check MyBasis Basis2 is expressed wrt CanonicalBasis
    @test b2.i.c0 == [0.0,1.0,0.0]
    @test b2.j.c0 == [0.0,0.0,1.0]
    @test b2.k.c0 == [1.0,0.0,0.0]
    @test b2.omega.c0 == [3.0,2.0,0.0]
    @test b2.alpha.c0 == [1.0,0.0,-5.0]

    # Check MyPoint Point 1 is expressed wrt CanonicalReferenceFrame
    @test p1.position.c0 == [1.0,4.0,2.0]
    @test p1.velocity.c0 == [2.0,0.0,13.0]
    @test p1.acceleration.c0 == [35.0,1.0,-6.0]

    # Check MyPoint Point 2 is expressed wrt CanonicalReferenceFrame
    @test p2.position.c0 == [1.0,4.0,-2.0]
    @test p2.velocity.c0 == [4.0,17.0,6.0]
    @test p2.acceleration.c0 == [8.0,-11.0,84.0]

    # Check Origin 1 expressed wrt rf1 is zero
    posOrigin1_1, velOrigin1_1, accOrigin1_1 = pos_vel_acc_inRF(origin1, rf1)
    @test posOrigin1_1 == [0.0,0.0,0.0]
    @test velOrigin1_1 == [0.0,0.0,0.0]
    @test accOrigin1_1 == [0.0,0.0,0.0]

    # Check Origin 2 expressed wrt rf2 is zero
    posOrigin2_2, velOrigin2_2, accOrigin2_2 = pos_vel_acc_inRF(origin2, rf2)
    @test posOrigin2_2 == [0.0,0.0,0.0]
    @test velOrigin2_2 == [0.0,0.0,0.0]
    @test accOrigin2_2 == [0.0,0.0,0.0]

    # Express CanonicalOriginPoint wrt rf2 
    origin_posRF2, origin_velRF2, origin_accRF2 = pos_vel_acc_inRF(CanonicalOriginPoint, rf2)
    @test origin_posRF2 == [-3.0,3.0,0.0]
    @test origin_velRF2 == [9.0,4.0,-8.0]
    @test origin_accRF2 == [0.0,-20.0,16.0]

    # Express Point 1 wrt rf2 using pos_vel_acc_inRF
    P1_posRF2, P1_velRF2, P1_accRF2 = pos_vel_acc_inRF(p1, rf2)
    @test P1_posRF2 == [1.0,5.0,1.0]
    @test P1_velRF2 == [15.0,7.0,-10.0]
    @test P1_accRF2 == [56.0,-48.0,-1.0]

    O1_posRF2, O1_velRF2, O1_accRF2 = pos_vel_acc_inRF(origin1, rf2)
    @test O1_posRF2 == [0.0,4.0,0.0]
    @test O1_velRF2 == [12.0,0.0,-10.0]
    @test O1_accRF2 == [4.0,-36.0,-1.0]

    # Find distance, relative velocity, and relative acceleration between particle2 and particle1 
    dist_P2P1_rf2 = norm(P1_posRF2 - posP2_22)
    dist_P2P1 = norm(p1.position.c0 - p2.position.c0)
    @test dist_P2P1_rf2 == 4.0
    @test dist_P2P1 == dist_P2P1_rf2
    relVel_P2P1_rf2 = velP2_22 - P1_velRF2
    relVel_P2P1_rf0 = p2.velocity.c0 - p1.velocity.c0
    @test relVel_P2P1_rf2 == [5.0,-7.0,10.0]
    @test relVel_P2P1_rf0 == [2.0,17.0,-7.0]
    #@test norm(relVel_P2P1_rf0) == norm(relVel_P2P1_rf2) # Check that the relative velocity is the same magnitude in both RF
    relAcc_P2P1_rf2 = accP2_22 - P1_accRF2
    relAcc_P2P1_rf0 = p2.acceleration.c0 - p1.acceleration.c0
    @test relAcc_P2P1_rf2 == [-58.0,48.0,1.0]
    @test relAcc_P2P1_rf0 == [-27.0,-12.0,90.0]
    #@test norm(relAcc_P2P1_rf0) == norm(relAcc_P2P1_rf2) # Check that the relative acceleration is the same magnitude in both RF
    
end


# Run the tests
@testset "LinearAlgebraTypes tests" begin 
    test_MyVector()
    test_MyBasis()
    test_MyPoint()
    test_MyReferenceFrame()
    test_MyParticle()
    test_componentsInBasis()
    test_pos_vel_acc_inRF()
    test_LinearAlgebraTypes()
    test_MyVector_operations()
    test_rotation_matrix()
    test_pure_rotation_MyBasis()
end