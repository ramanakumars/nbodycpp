/**
 * @file hermite.cpp
 * @brief Implementation of Hermite 4th order predictor-corrector integrator
 *
 * The Hermite method is widely used in astrophysical N-body simulations for its
 * excellent energy conservation and ability to handle close encounters.
 *
 * Key advantages:
 * - 4th order accuracy with only ~2 force evaluations per step
 * - Superior energy conservation compared to symplectic integrators
 * - Better handling of close encounters
 * - Can use 2-3x larger timesteps than Yoshida for same accuracy
 *
 * Reference: Makino & Aarseth (1992), PASJ, 44, 141
 */

#include "hermite.h"
#include "barneshut.h"

void getAccelerationAndJerk(std::vector<std::shared_ptr<Particle>> &particles,
                           QuadTree<Particle> *tree)
{
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++)
    {
        particles[i]->acceleration.zero();
        particles[i]->jerk.zero();
        BarnesHutForceAndJerk(particles[i].get(), tree, 0.05);
    }
}

/**
 * @brief One timestep of Hermite 4th order integration
 *
 * @details Three-stage predictor-corrector scheme:
 *
 * **Stage 1 - Predictor:**
 * Extrapolate to t+dt using Taylor series with jerk:
 * - x_p = x + v·dt + ½a·dt² + ⅙j·dt³
 * - v_p = v + a·dt + ½j·dt²
 *
 * **Stage 2 - Evaluator:**
 * Calculate forces at predicted positions
 *
 * **Stage 3 - Corrector:**
 * Update using average of old and new derivatives:
 * - v_new = v + ½(a₀+a₁)·dt + 1/12(j₀-j₁)·dt²
 * - x_new = x + ½(v₀+v₁)·dt + 1/12(a₀-a₁)·dt²
 *
 * @param particles Particle array to integrate
 * @param tree QuadTree for force calculation
 * @param dt Timestep size
 *
 * @note Requires particles to have jerk initialized from previous step
 * @note Thread-safe with OpenMP parallelization
 * @note Cost: ~2 force evaluations (predictor + corrector)
 */
void hermiteStep(std::vector<std::shared_ptr<Particle>> &particles,
                QuadTree<Particle> *tree, double dt)
{
    // PREDICTOR: Predict positions and velocities at t + dt
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++)
    {
        auto &p = particles[i];

        // 4th order predictor
        // x_pred = x + v*dt + (1/2)*a*dt² + (1/6)*jerk*dt³
        // v_pred = v + a*dt + (1/2)*jerk*dt²

        p->position_pred = p->position +
                          p->velocity * dt +
                          p->acceleration * (0.5 * dt * dt) +
                          p->jerk * (dt * dt * dt / 6.0);

        p->velocity_pred = p->velocity +
                          p->acceleration * dt +
                          p->jerk * (0.5 * dt * dt);
    }

    // EVALUATOR: Calculate forces and jerks at predicted positions
    // Temporarily swap positions for force calculation
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++)
    {
        auto &p = particles[i];
        std::swap(p->position, p->position_pred);
        std::swap(p->velocity, p->velocity_pred);
    }

    // Store old acceleration and jerk
    std::vector<vector2D> old_acc(particles.size());
    std::vector<vector2D> old_jerk(particles.size());

#pragma omp parallel for schedule(static, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++)
    {
        old_acc[i] = particles[i]->acceleration;
        old_jerk[i] = particles[i]->jerk;
    }

    // Calculate new accelerations and jerks at predicted positions
    getAccelerationAndJerk(particles, tree);

    // CORRECTOR: Correct positions and velocities using improved estimates
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++)
    {
        auto &p = particles[i];

        vector2D a0 = old_acc[i];
        vector2D a1 = p->acceleration;
        vector2D jerk0 = old_jerk[i];
        vector2D jerk1 = p->jerk;

        // Swap back to get original positions/velocities
        std::swap(p->position, p->position_pred);
        std::swap(p->velocity, p->velocity_pred);

        // 4th order corrector using both old and new derivatives
        // v_new = v_old + (a0 + a1)*dt/2 + (jerk0 - jerk1)*dt²/12
        // x_new = x_old + (v_old + v_new)*dt/2 + (a0 - a1)*dt²/12

        vector2D v_old = p->velocity;

        p->velocity = v_old +
                     (a0 + a1) * (0.5 * dt) +
                     (jerk0 - jerk1) * (dt * dt / 12.0);

        p->position = p->position +
                     (v_old + p->velocity) * (0.5 * dt) +
                     (a0 - a1) * (dt * dt / 12.0);
    }
}
