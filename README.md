![AstrodynamicsEdu Logo](docs/assets/LogoAstrodynamicsEdu.png | width=100)

# AstrodynamicsEdu
_The AstrodynamicsEdu.jl is an educational and user-friendly astrodynamics package implemented in Julia. It provides a comprehensive set of tools and algorithms for performing various calculations and simulations related to orbital mechanics, as well as, the visulization of orbits in both 2D and 3D. With AstrodynamicsEdu.jl, undergraduate students, researchers, and enthusiasts alike can gain a deeper understanding of Orbital Dynamics._

[![Build Status](https://github.com/AliciaSBa/Astrodynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AliciaSBa/Astrodynamics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://aliciasba.github.io/AstrodynamicsEdu.jl/)

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
        <li><a href="#features">Features</a></li>
        <li><a href="#built-with">Built With</a></li>
      </ul>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About The Project

Nowadays, computational tools and libraries have become essential in solving Astrodynamics problems. However, there is a significant gap in the availability of educational resources and software libraries in Julia that specifically address Orbital Dynamics. The aim of this library is to develop an educational Astrodynamics library in Julia which is accessible to students with a basic knowledge of programming and a wish to learn Orbital Dynamics in a dynamic and visual way.


<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Features

- **Solving Problems of the Ideal Two-Body Problem:** This package implements numerical algorithms to solve problems related to the Ideal Two-Body Problem. This involves modeling the motion of a spacecraft under the influence of gravitational forces from a central body. This allows to compute and analyze the spacecraft’s position and velocity at any given time, facilitating the study of various orbital phenomena and the prediction of future orbits.

- **Change of Reference Frame:** AstrodynamicsEdu.jl supports the transformation of orbital elements and state vectors between different reference frames, such as inertial frames, rotating frames, and non-inertial frames. Which allows to analyze orbits from different perspectives and understand the impact of reference frame choice on orbital parameters and dynamics.

- **Orbital Maneuvers and Delta-v Computation:** Functions are provided to compute common orbital maneuvers, such as Hohmann transfers and inclination changes. Students will be able to input the initial and target orbits, and the code will calculate the required delta-v (change in velocity) for executing the maneuvers. This aids in understanding the principles behind orbital maneuvers and helps to gain practical experience in mission planning.

- **Orbital Perturbations:** The package incorporates models for common orbital perturbations such as atmospheric drag and Earth's oblateness. The user has the option to include these perturbations in their simulations and observe the effects on spacecraft trajectories.

- **Orbit Visualization and Plotting:** The code includes functionality to plot and visualize orbits in two or three dimensions. The user will be able to visualize spacecraft trajectories, observe the effects of different orbital parameters, and analyze orbital elements such as semimajor axis, eccentricity, inclination, and argument of periapsis. These visualizations enhance the understanding of orbital dynamics concepts and aid in problem-solving exercises.

- **Educational Resources and Examples:** This package also provides educational resources, including sample problems and examples, to help students reinforce their understanding of orbital dynamics concepts. These resources will cover a range of topics such as Kepler’s laws, orbital elements, vis-viva equation, and orbital transfers.


<p align="right">(<a href="#readme-top">back to top</a>)</p>

### Built With
Here are listed the major libraries that have been used to bootstrap this project.
* LinearAlegbra.jl
* Plots.jl
* Roots.jl
* DifferentialEquations.jl
* and more

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->
## Getting Started

### Prerequisites

To use this library, first you will need to have:
* **Julia 1.8.3**. To install this version of the Julia software follow the instructions explained in the following link: [https://julialang.org/downloads/](https://julialang.org/downloads/)


### Installation

This package can be installed using:
  ```sh
  julia> Pkg.update()
  julia> Pkg.add("AstrodynamicsEdu")
  ```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

Some quick examples can be seen below!
  ```sh
 # Installation
 import Pkg
 Pkg.add("AstrodynamicsEdu") # or julia> ]install Astrodynamics
 
  # Loading
  using AstrodynamicsEdu
  
  # Solve the I2BP 
  ```

_For more examples, please see the [Package Documentation](https://example.com)_

<p align="right">(<a href="#readme-top">back to top</a>)</p>
