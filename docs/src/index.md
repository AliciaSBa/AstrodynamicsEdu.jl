# AstrodynamicsEdu.jl Documentation

Welcome to the documentation for AstrodynamicsEdu, an educational package designed to provide 
a comprehensive understanding of orbital dynamics. This documentation serves as a guide for 
students, educators, and enthusiasts interested in exploring the fascinating world of astrodynamics.

AstrodynamicsEdu offers a wide range of features to facilitate learning and problem-solving in the 
field of orbital dynamics. Let's take a closer look at the key features offered by this package:

1. Solving Problems of the Ideal Two-Body Problem.
2. Change of Reference Frame: AstrodynamicsEdu supports the transformation of orbital elements and state vectors between different reference frames, including inertial frames, rotating frames, and non-inertial frames.
3. Orbital Maneuvers and Delta-v Computation.
4. Orbital Perturbations: AstrodynamicsEdu incorporates models for common perturbations, including atmospheric drag and Earth's oblateness.
5. Orbit Visualization and Plotting.
6. Educational Resources and Examples: To support the learning process, AstrodynamicsEdu provides a collection of educational resources, sample problems, and examples.

We hope that AstrodynamicsEdu and this documentation empower you on your journey to mastering orbital 
dynamics.

## Getting Started: Installation And First Steps
As this package has not been registered yet, to run this package you need to write the following 
commands in your Julia REPL terminal to install the package:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.add(url="https://github.com/AliciaSBa/AstrodynamicsEdu.jl.git")
```

Once it has been registered, it will only be necessary to run the following command to install it:
```julia
Pkg.add("AstrodynamicsEdu")
```

To load the package, use the command:

```julia
using AstrodynamicsEdu
```

## Package Structure
The AstrodynamicsEdu package is composed by six different submodules, each designed to provide a 
different set of capabilities related with Orbital Dynamics. The different functionalities and 
details of each specific module can be found on its respective documentation page in the Modules section.
```@contents
Pages = [
"modules/linearAlgebraTypes.md",
"modules/idealTwoBodyProblem.md"
]
Depth = 2
```

## Examples
To help you kick-start using the AstrodynamicsEdu package, a set of extensive orbital dynamics
problems are provided as examples. These problems have been taken from the course of Space Vehicles 
and Orbital Dynamics from University Carlos III of Madrid (uc3m). You can find the example 
problems in the Tutorials section.

```@contents
Pages = [
    "tutorials/example1.md",
]
Depth = 2
```