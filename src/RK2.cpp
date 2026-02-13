/**
 * @file RK2.cpp
 * @brief Implementation of second-order Runge-Kutta integrator
 */

#include "RK2.h"
#include "barneshut.h"

/**
 * @brief RK2 integration step (midpoint method)
 *
 * @details Steps:
 * 1. Calculate current acceleration
 * 2. Predictor: advance to midpoint using current velocity and acceleration
 * 3. Evaluate acceleration at predicted position
 * 4. Corrector: update velocity using average of initial and midpoint accelerations
 *
 * @param particles Particle vector
 * @param tree QuadTree for force calculation
 * @param dt Timestep
 */
void RK2step(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree,
             double dt) {
    getAcceleration(particles, tree);

#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++) {
        auto &p = particles[i];
        Particle temp = *p;

        // First half RK2 step
        temp.position += p->velocity * dt + p->acceleration * (0.5 * dt * dt);

        // Get intermediate acceleration
        temp.acceleration.zero();
        BarnesHutForceAndJerk(&temp, tree, 0.05);

        // Second half RK2 step
        p->position = temp.position;
        p->velocity += (temp.acceleration + p->acceleration) * (0.5 * dt);
        p->acceleration = temp.acceleration;
    }
}
