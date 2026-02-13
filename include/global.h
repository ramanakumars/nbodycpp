/**
 * @file global.h
 * @brief Global constants, includes, and type definitions for N-body simulation
 */

#pragma once

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <memory>
#include <iostream>
#include <algorithm>
#include <unordered_set>
#include <omp.h>

/// @brief OpenMP chunk size for parallel loops
#define CHUNK_SIZE 500

/// @brief Gravitational constant (in simulation units)
const double GRAV_G = 1;

/// @brief Reference mass for Barnes-Hut theta scaling
const double MASS_REF = 0.1;

/// @brief Exponent for theta scaling: theta_eff = theta * (M_ref/M_cell)^alpha
const double ALPHA = 0.5;

/**
 * @enum transport_type
 * @brief Available integrator types
 */
enum transport_type {
    RK2,     ///< 2nd order Runge-Kutta (midpoint method)
    YOSHIDA, ///< 4th order symplectic Yoshida integrator
    HERMITE  ///< 4th order Hermite predictor-corrector
};

/// @brief Global integrator selection
extern transport_type TRANSPORT_TYPE;
