# Linear Algebra Types

In order to be able to perform operations more at ease and to store information easily, different linear algebra types were created. Their main feature is that through different constructors the user can pass the components they desire and these will always be stored in the Canonical Reference Frame or Basis. Having all vectors, points, vector bases, and reference frame components stored in only one basis/reference frame allows to speed the computation and to reduce the amount of possible errors.

The diagram below shows the five different types defined in the linear algebra module. Below each type, it can be seen the objects they store and the type of each of them (indicated as \textit{fieldName::Type}). To access the field value inside the object it is only necessary to write the name of the variable followed by the fieldName. For example, if 'B1' is MyBasis object and I want to access its angular velocity, I would just have to type 'B1.omega'.

![Alt text for screen readers](https://github.com/AliciaSBa/AstrodynamicsEdu.jl/edit/main/docs/src/modules/LinearAlgebraTypes.png "Linear Algebra Types Diagram")

All these linear algebra types are related among each other, that is why it was so important to get the chain of dependency between them correctly. This is where **constructors** come in. The Julia language allows to have more than one way to build a type/struct. Thanks to this feature, the user can have a vector v₁ expressed in vector basis B₁, and if both v₁ and B₁ (which must be a MyBasis, this means, it must be expressed with respect to the Canonical Basis) are passed, v₁ can be turned into a MyVector. This means simply that v₁ is stored internally with the coordinates with respect to the Canonical Basis. So, now, we will examine which are the different available constructors for each of the types.

1. **MyVector**: stores a vector of length three in the canonical right-handed vector basis B₀. These are its possible constructors:
   - **MyVector(c₀)**, where c₀ is a Vector{Float64} of length = 3 already expressed in the Canonical Basis.
   - **MyVector(c₁, B₁)**, where c₁ is a Vector{Float64} expressed in vector basis B₁, which is a MyBasis object.
   - **MyVector(x₀, y₀, z₀)**, where x₀, y₀, and z₀ are the vector components (numbers) already expressed in the Canonical Basis.
   - **MyVector(x₁, y₁, z₁, B₁)**, where x₁, y₁, and z₁ are the vector components (numbers) expressed in vector basis B₁, which is a MyBasis object.
   - **MyVector()**, which simply builds a MyVector object full of zeros.

2. **MyBasis**: defines a vector basis (i.e., B₁) by storing its unit vectors, angular velocity, and angular acceleration with respect to the canonical basis B₀. These are the possible constructors:
   - **MyBasis(i₁₀, j₁₀, k₁₀, omega₁₀, alpha₁₀)**, where all are MyVector objects. i₁₀, j₁₀, k₁₀ are the unit vectors of vector basis B₁, and omega₁₀ and alpha₁₀ are the angular velocity and angular acceleration, respectively.
   - **MyBasis(i₂₁, j₂₁, k₂₁, omega₂₁, alpha₂₁, B₁)**, where all are Vector{Float64} which express the vector basis B₂, with respect to the vector basis B₁, which is a MyBasis object.
   - **MyBasis()**, which builds the Canonical Basis.

3. **MyPoint**: stores the position, velocity, and acceleration vectors that define a point with respect to the inertial reference frame S₀. These are the possible constructors:
   - **MyPoint(pos, veloc, accel)**, where pos, veloc, and accel are the position, velocity, and acceleration vectors already MyVector objects as they are expressed with respect to the Canonical Reference Frame.
   - **MyPoint(pos₁, veloc₁, accel₁, RF₁)**, where pos₁, veloc₁, and accel₁ are Vector{Float64} expressed with respect to the reference frame RF₁, which is a MyReferenceFrame object.
   - **MyPoint()**, which builds the Canonical Origin Point. This is just an empty point in the Canonical Reference Frame.

4. **MyReferenceFrame**: defined by a vector basis (i.e., B₁) which is a MyBasis object and a point (i.e., P₁) which is a MyPoint object. These are the possible constructors:
   - **MyReferenceFrame(basis, origin)**, where basis is a MyBasis object and origin is a MyPoint object.
   - **MyReferenceFrame()**, which returns the inertial reference frame, whose basis is the Canonical Basis and origin point is the Canonical Origin Point.

5. **MyParticle**: object that defines the point particle by storing the point, which is a MyPoint object, and its associated mass. This is the constructor:
   - **MyParticle(point, mass)**, where point is a MyPoint object and mass is just a number.



Some useful objects have been also predefined, so that the user can simply call them every time they are needed, allowing therefore to save time. These **predefined objects** are the canonical vector basis $B_0$, the three canonical unitary vectors ($\mathbf{i_0},\mathbf{j_0},\mathbf{k_0}$), the canonical origin point $O_0$, and the canonical/inertial reference frame $S_0$. They can simply be used by typing their correct name, which appear highlighted below:

- **CanonicalBasis** = MyBasis(MyVector([1.0, 0.0, 0.0]), MyVector([0.0, 1.0, 0.0]), MyVector([0.0, 0.0, 1.0]), MyVector([0.0, 0.0, 0.0]), MyVector([0.0, 0.0, 0.0]))
- **i0** = MyVector([1.0, 0.0, 0.0])
- **j0** = MyVector([0.0, 1.0, 0.0])
- **k0** = MyVector([0.0, 0.0, 1.0])
- **CanonicalOriginPoint** = MyPoint(MyVector([0.0, 0.0, 0.0]), MyVector([0.0, 0.0, 0.0]), MyVector([0.0, 0.0, 0.0]))
- **CanonicalReferenceFrame** = MyReferenceFrame(MyBasis(MyVector([1.0, 0.0, 0.0]), MyVector([0.0, 1.0, 0.0]), MyVector([0.0, 0.0, 1.0]), MyVector([0.0, 0.0, 0.0]), MyVector([0.0, 0.0, 0.0])), MyPoint(MyVector([0.0, 0.0, 0.0]), MyVector([0.0, 0.0, 0.0]), MyVector([0.0, 0.0, 0.0])))


Both recipes and functions have been implemented to allow plotting MyBasis and MyReferenceFrame objects in 3D. The only difference in their usage is that for the recipes it is required to include again the Plots.jl package, by writing before:_using Plots_. The proper way to call this functions or recipes is addressed in the table below:

| Function call               | Recipe call  | Functionality                                                                                  |
|-----------------------------|--------------|------------------------------------------------------------------------------------------------|
| plot_MyBasis(B1)            | plot(B1)       | Plot the MyBasis object B1 in 3D using either the function or the recipe                       |
| plot_MyReferenceFrame(RF1)  | plot(RF1)      | Plot the MyReferenceFrame object S1 in 3D using either the function or the recipe             |


In addition, some useful functions have also been created to assist in changing of
reference frame, rotating vector basis, and more. 

| Function Name            | Inputs                                                | Outputs                                      | Functionality                                                                                                |
|--------------------------|-------------------------------------------------------|----------------------------------------------|--------------------------------------------------------------------------------------------------------------|
| componentsInBasis        | v::MyVector<br>B1::MyBasis                            | v1::Vector                                   | Project a MyVector object onto a basis B₁                                                                      |
| pos_vel_acc_inRF         | point::MyPoint<br>RF1::MyReferenceFrame               | pos1::Vector<br>veloc1::Vector<br>accel1::Vector | Obtain the position, velocity, and acceleration vectors of a MyPoint object in a different reference frame S₁ |
| rotation_matrix          | basis::MyBasis                                        | R::Matrix                                    | Obtain the Rotation Matrix [₀R₁] of a MyBasis object B₁                                                        |
| pure_rotation_MyBasis_x  | basis::MyBasis<br>theta::Number                       | rotB::MyBasis                               | Obtain the new MyBasis object after a pure rotation about the x-axis                                         |
| pure_rotation_MyBasis_y  | basis::MyBasis<br>theta::Number                       | rotB::MyBasis                               | Obtain the new MyBasis object after a pure rotation about the y-axis                                         |
| pure_rotation_MyBasis_z  | basis::MyBasis<br>theta::Number                       | rotB::MyBasis                               | Obtain the new MyBasis object after a pure rotation about the z-axis                                         |
