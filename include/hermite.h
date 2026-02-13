/**
 * @file hermite.h
 * @brief Hermite 4th order predictor-corrector integrator for N-body simulations
 *
 * This file implements the Hermite integration scheme, which is widely used
 * in astrophysical N-body simulations due to its superior energy conservation
 * and handling of close encounters compared to symplectic integrators.
 *
 * The Hermite method uses both force (acceleration) and force derivative (jerk)
 * to achieve 4th order accuracy with only 2 force evaluations per timestep.
 */

#pragma once

#include "global.h"
#include "particle.h"
#include "quadtree.h"

/**
 * @brief Performs one timestep using the Hermite 4th order predictor-corrector method
 *
 * @details The Hermite integrator consists of three stages:
 *
 * 1. **Predictor**: Predicts positions and velocities at t+dt using Taylor expansion:
 *    - x_pred = x + v*dt + (1/2)*a*dt² + (1/6)*jerk*dt³
 *    - v_pred = v + a*dt + (1/2)*jerk*dt²
 *
 * 2. **Evaluator**: Calculates forces and jerks at predicted positions
 *
 * 3. **Corrector**: Corrects positions/velocities using average of old and new derivatives:
 *    - v_new = v_old + (a0 + a1)*dt/2 + (jerk0 - jerk1)*dt²/12
 *    - x_new = x_old + (v_old + v_new)*dt/2 + (a0 - a1)*dt²/12
 *
 * @param particles Vector of particles to integrate
 * @param tree QuadTree for Barnes-Hut force calculation
 * @param dt Timestep size
 *
 * @note Requires particles to have jerk field initialized
 * @note More accurate than Yoshida-4 for same timestep, with ~2x force evaluations
 * @note Non-symplectic but excellent energy conservation in practice
 */
void hermiteStep(std::vector<std::shared_ptr<Particle>> &particles,
                 QuadTree<Particle> *tree, double dt);

/**
 * @brief Calculates both acceleration and jerk (time derivative of acceleration) for all particles
 *
 * @details Uses Barnes-Hut algorithm for O(N log N) complexity.
 * The jerk is calculated as the time derivative of gravitational acceleration:
 *
 * For particle pair: da/dt = -G*m * [v/r³ - 3*(r·v)*r/r⁵]
 *
 * where r is position difference and v is velocity difference.
 *
 * @param particles Vector of particles
 * @param tree QuadTree for hierarchical force calculation
 *
 * @note This function zeros acceleration and jerk before calculation
 * @note Parallelized with OpenMP
 */
void getAccelerationAndJerk(std::vector<std::shared_ptr<Particle>> &particles,
                           QuadTree<Particle> *tree);

