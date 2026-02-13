/**
 * @file yoshida.h
 * @brief Yoshida 4th order symplectic integrator
 *
 * Implements the Yoshida symplectic integrator, a composition method that
 * achieves 4th order accuracy while preserving the symplectic structure
 * of Hamiltonian systems.
 *
 * This results in excellent long-term energy conservation, making it ideal
 * for orbital dynamics and secular evolution studies.
 *
 * Reference: Yoshida (1990), Physics Letters A, 150, 262
 */

#pragma once

#include "global.h"
#include "barneshut.h"

/// @name Yoshida coefficients
/// @{
const double w0 = -pow(2.0, 1.0 / 3.0) / (2.0 - pow(2.0, 1.0 / 3.0)); ///< Weight 0
const double w1 = 1.0 / (2.0 - pow(2.0, 1.0 / 3.0));                  ///< Weight 1

const double c1 = w1 / 2.0;          ///< Drift coefficient 1
const double c2 = (w0 + w1) / 2.0;   ///< Drift coefficient 2
const double c3 = c2;                ///< Drift coefficient 3
const double c4 = c1;                ///< Drift coefficient 4

const double d1 = w1;                ///< Kick coefficient 1
const double d2 = w0;                ///< Kick coefficient 2
const double d3 = w1;                ///< Kick coefficient 3
/// @}

/**
 * @brief Drift step: update positions using current velocities
 *
 * @details x_new = x_old + v * dt
 *
 * @param particles Particle vector
 * @param dt Timestep (may be fractional for substeps)
 */
void drift(std::vector<Particle> &, double);

/**
 * @brief Kick step: update velocities using current accelerations
 *
 * @details v_new = v_old + a * dt
 *
 * @param particles Particle vector
 * @param dt Timestep (may be fractional for substeps)
 */
void kick(std::vector<Particle> &, double);

/**
 * @brief Perform one timestep using Yoshida 4th order symplectic method
 *
 * @details Composition of drift-kick-drift stages with special coefficients:
 * - Drift(c1*dt) → Kick(d1*dt) → Drift(c2*dt) → ...
 * - Total of 4 drift and 3 kick substeps
 * - Coefficients chosen to cancel error terms up to 4th order
 *
 * Properties:
 * - Symplectic (preserves phase space volume)
 * - Time-reversible
 * - 4th order accuracy: error ~ O(dt⁵)
 * - Excellent long-term energy conservation
 *
 * @param particles Vector of particles to integrate
 * @param tree QuadTree for Barnes-Hut force calculation
 * @param dt Timestep size
 *
 * @note Requires 4 force evaluations per step
 * @note Best for long-term orbital integration
 * @note May struggle with close encounters
 */
void yoshidaStep(std::vector<std::shared_ptr<Particle>> &, QuadTree<Particle> *, double);