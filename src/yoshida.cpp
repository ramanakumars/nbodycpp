/**
 * @file yoshida.cpp
 * @brief Implementation of Yoshida 4th order symplectic integrator
 */

#include "yoshida.h"

/**
 * @brief Drift: advance positions by velocity
 */
void drift(std::vector<std::shared_ptr<Particle>> &particles, double dt) {
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++) {
        auto &p = particles[i];
        p->position += p->velocity * dt;
    }
}

/**
 * @brief Kick: advance velocities by acceleration
 */
void kick(std::vector<std::shared_ptr<Particle>> &particles, double dt) {
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++) {
        auto &p = particles[i];
        p->velocity += p->acceleration * dt;
    }
}

/**
 * @brief Yoshida 4th order integration step
 *
 * @details Performs composition of drift and kick operations with
 * specially chosen coefficients for 4th order accuracy while maintaining
 * symplectic properties.
 */
void yoshidaStep(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree,
                 double dt) {
    // First stage
    drift(particles, c1 * dt);
    getAcceleration(particles, tree);
    kick(particles, d1 * dt);

    // Second stage
    drift(particles, c2 * dt);
    getAcceleration(particles, tree);
    kick(particles, d2 * dt);

    // Third stage
    drift(particles, c3 * dt);
    getAcceleration(particles, tree);
    kick(particles, d3 * dt);

    drift(particles, c4 * dt);
}
