![AstrodynamicsEdu Logo](https://www.dropbox.com/s/avv0cej5t5ljp8e/LogoAstrodynamicsEdu.png?dl=0)

# AstrodynamicsEdu
_Astrodynamics calculations and orbits visualization to solve Orbital Dynamics problems._

[![Build Status](https://github.com/AliciaSBa/Astrodynamics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AliciaSBa/Astrodynamics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://aliciasba.github.io/AstrodynamicsEdu.jl/)

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
      <ul>
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
