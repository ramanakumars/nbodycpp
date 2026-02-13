/**
 * @file RK2.h
 * @brief Second-order Runge-Kutta (RK2) integrator
 *
 * Also known as the midpoint method or modified Euler method.
 * Provides 2nd order accuracy with 2 force evaluations per step.
 */

#pragma once

#include "barneshut.h"
#include "global.h"

/**
 * @brief Performs one timestep using 2nd order Runge-Kutta method
 *
 * @details The RK2 method is a predictor-corrector scheme:
 * 1. Predict position at midpoint: x_mid = x + v*dt + 0.5*a*dt²
 * 2. Evaluate acceleration at midpoint
 * 3. Update velocity using average: v_new = v + 0.5*(a_0 + a_mid)*dt
 * 4. Update position to x_mid
 *
 * This provides better accuracy than Euler's method but is not symplectic,
 * so it may have energy drift over long integrations.
 *
 * @param particles Vector of particles to integrate
 * @param tree QuadTree for Barnes-Hut force calculation
 * @param dt Timestep size
 *
 * @note Non-symplectic (not ideal for long-term orbit integration)
 * @note Cost: 2 force evaluations per step
 * @note Better suited for short-term high-accuracy calculations
 */
void RK2step(std::vector<std::shared_ptr<Particle>> &, QuadTree<Particle> *, double);
