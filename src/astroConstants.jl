# AstroConstants file

export G, k, sigma, c, 
        mu_Sun, mu_Mercury, mu_Venus, mu_Earth, mu_Moon, mu_Mars, mu_Jupiter, mu_Saturn, mu_Uranus, mu_Neptune, 
        R_Sun, R_Mercury, R_Venus, R_Earth, R_Moon, R_Mars, R_Jupiter, R_Saturn, R_Uranus, R_Neptune, 
        e_Mercury, e_Venus, e_Earth, e_Moon, e_Mars, e_Jupiter, e_Saturn, e_Uranus, e_Neptune, 
        a_Mercury, a_Venus, a_Earth, a_Moon, a_Mars, a_Jupiter, a_Saturn, a_Uranus, a_Neptune, 
        omega_Mercury, omega_Venus, omega_Earth, omega_Moon, omega_Mars, omega_Jupiter, omega_Saturn, omega_Uranus, omega_Neptune,
        J2_Mercury, J2_Venus, J2_Earth, J2_Moon, J2_Mars, J2_Jupiter, J2_Saturn, J2_Uranus, J2_Neptune,
        M_Mercury, M_Venus, M_Earth, M_Moon, M_Mars, M_Jupiter, M_Saturn, M_Uranus, M_Neptune

# PHYSICAL CONSTANTS

# Gravitational constant, G [m^3/(kg*s^2)]
const G = 6.67408e-11

# Boltzmann constant, k [m^2*kg/(s^2*K)]
const k = 1.38064852e-23

# Stefan-Boltzmann constant, sigma [W/(m^2*K^4)]
const sigma = 5.670367e-8

# Speed of light in vacuum, c [m/s]
const c = 299792458.0


# CELESTIAL BODY DATA 
# (from NASA's planetary fact sheet - Date: 11-02-2023)

# Gravitational parameter, mu [km^3/s^2]
const mu_Sun = 1.32712440018e11
const mu_Mercury = 2.2032e4
const mu_Venus = 3.24859e5
const mu_Earth = 398600.4418
const mu_Moon = 4902.800066
const mu_Mars = 4.282837e4
const mu_Jupiter = 1.26687e8
const mu_Saturn = 3.7931187e7
const mu_Uranus = 5.79395e6
const mu_Neptune = 6.836510e6

# Mean body radius, R [km]
const R_Sun = 6.957e5
const R_Mercury = 2.4397e3
const R_Venus = 6.0518e3
const R_Earth = 6.371e3
const R_Moon = 1.7374e3
const R_Mars = 3.3895e3
const R_Jupiter = 6.9911e4
const R_Saturn = 5.8232e4
const R_Uranus = 2.5362e4
const R_Neptune = 2.4622e4

# Eccentricity, e [-]
const e_Mercury = 0.2056
const e_Venus = 0.0068
const e_Earth = 0.0167
const e_Moon = 0.0549
const e_Mars = 0.0935
const e_Jupiter = 0.0487
const e_Saturn = 0.0520
const e_Uranus = 0.0469
const e_Neptune = 0.0097

# Semi-major axis, a [km]
const a_Mercury = 5.790905e7
const a_Venus = 1.08209e8
const a_Earth = 1.49598e8
const a_Moon = 3.8475e5
const a_Mars = 2.2796e8
const a_Jupiter = 7.7848e8
const a_Saturn = 1.4320e9
const a_Uranus = 2.8670e9
const a_Neptune = 4.5150e9

# Sidereal rotation rate of the body about its axis, omega [rad/s]
const omega_Mercury = 1.240014e-6
const omega_Venus = -0.29924e-6
const omega_Earth = 7.292115e-5
const omega_Moon = 2.6617e-6
const omega_Mars = 7.088218e-5
const omega_Jupiter = 1.7585e-4
const omega_Saturn = 1.63785e-4
const omega_Uranus = -1.01237e-4
const omega_Neptune = 1.08338e-4

# Oblateness parameter of the body's gravitational field, J2 [-]
const J2_Mercury = 5.03e-5
const J2_Venus = 4.458e-6
const J2_Earth = 1.08263e-3
const J2_Moon = 2.027e-4
const J2_Mars = 1.96045e-3
const J2_Jupiter = 1.4736e-2
const J2_Saturn = 1.6298e-2
const J2_Uranus = 3.34343e-3
const J2_Neptune = 3.411e-3

# Mass of the body, M [kg]
const M_Sun = 1.988544e30
const M_Mercury = 3.302e23
const M_Venus = 4.8685e24
const M_Earth = 5.97219e24
const M_Moon = 7.34767309e22
const M_Mars = 6.4171e23
const M_Jupiter = 1.89813e27
const M_Saturn = 5.68319e26
const M_Uranus = 8.68103e25 
const M_Neptune = 1.0241e26

