/**
 * @file interactions.cpp
 * @brief Implementation of particle interactions and collision handling
 */

#include "RK2.h"
#include "hermite.h"
#include "interactions.h"
#include "yoshida.h"

/// @brief Selected integrator (default: Hermite 4th order)
transport_type TRANSPORT_TYPE = transport_type::HERMITE;

/**
 * @brief Main particle update: integrate and handle collisions
 */
void updateParticles(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree,
                     double dt) {
    transportStep(particles, tree, dt);

    checkCollisions(particles, tree, dt);

    double total_mass = 0, com_x = 0, com_y = 0;
#pragma omp parallel for reduction(+ : com_x, com_y, total_mass) schedule(static, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++) {
        if (tree->bounds.contains(particles[i]->position)) {
            com_x += particles[i]->position.x * particles[i]->mass;
            com_y += particles[i]->position.y * particles[i]->mass;
            total_mass += particles[i]->mass;
        }
    }

    com_x /= total_mass;
    com_y /= total_mass;

    vector2D new_com = {com_x, com_y};
#pragma omp parallel for schedule(static, CHUNK_SIZE)
    for (auto &particle : particles) {
        particle->position -= new_com;
    }
}

/**
 * @brief Dispatch to selected integrator
 */
void transportStep(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree,
                   double dt) {
    if (TRANSPORT_TYPE == YOSHIDA) {
        yoshidaStep(particles, tree, dt);
    } else if (TRANSPORT_TYPE == RK2) {
        RK2step(particles, tree, dt);
    } else if (TRANSPORT_TYPE == HERMITE) {
        hermiteStep(particles, tree, dt);
    } else {
        fprintf(stderr, "%d not a valid transport type\n", TRANSPORT_TYPE);
        exit(1);
    }
}

/**
 * @brief Calculate relative gravitational acceleration for collision prediction
 *
 * @details Computes acceleration of p1 relative to p2 in their mutual frame.
 * Used for trajectory prediction in continuous collision detection.
 *
 * @param p1 First particle
 * @param p2 Second particle
 * @return Relative acceleration vector
 */
vector2D calcMutualAcceleration(const Particle *p1, const Particle *p2) {
    vector2D diff = p1->position - p2->position;
    double dist = std::max(diff.norm(), p1->radius + p2->radius);
    double invDistCubed = 1.0 / (dist * dist * dist);

    // Relative acceleration: a1 - a2 = -G*(m1+m2)*r/|r|³
    double scale = -GRAV_G * (p1->mass + p2->mass) * invDistCubed;
    return diff * scale;
}

// Continuous collision detection with gravitational effects
// Uses quadratic trajectory approximation: pos(t) = p0 + v0*t + 0.5*a*t²
CollisionInfo predictCollision(const Particle *p1, const Particle *p2, double dt) {
    CollisionInfo info;

    vector2D relPos = p1->position - p2->position;
    vector2D relVel = p1->velocity - p2->velocity;

    // Calculate relative acceleration (includes gravity)
    vector2D relAcc = calcMutualAcceleration(p1, p2);

    double collisionRadius = p1->radius + p2->radius;

    // For very close particles, use simple distance check with current acceleration
    double currentDist = relPos.norm();
    if (currentDist < collisionRadius * 1.1) {
        info.willCollide = true;
        info.collisionTime = 0;
        info.minDistance = currentDist;
        return info;
    }

    // Trajectory with constant acceleration: r(t) = r0 + v0*t + 0.5*a*t²
    // We need to solve: |r(t)|² = R²
    // This gives a quartic equation, but we can find minimum distance using calculus

    // Find time of closest approach by minimizing |r(t)|²
    // d/dt[|r(t)|²] = 0
    // 2*r·(v + a*t) = 0
    // r·v + r·a*t + v·a*t + a·a*t² = 0
    // But this is approximate. For better accuracy, sample the trajectory.

    const int numSamples = 10;
    double minDist = currentDist;
    double minDistTime = 0;
    bool foundCollision = false;
    double collisionTime = dt;

    for (int i = 0; i <= numSamples; i++) {
        double t = (dt * i) / numSamples;

        // Position at time t with quadratic trajectory
        vector2D posAtT = relPos + relVel * t + relAcc * (0.5 * t * t);
        double distAtT = posAtT.norm();

        if (distAtT < minDist) {
            minDist = distAtT;
            minDistTime = t;
        }

        if (distAtT < collisionRadius && !foundCollision) {
            foundCollision = true;
            collisionTime = t;
        }
    }

    // Refine minimum distance estimate around the minimum
    if (minDistTime > 0 && minDistTime < dt) {
        double refineDt = dt / (2.0 * numSamples);
        for (int i = -2; i <= 2; i++) {
            double t = std::max(0.0, std::min(dt, minDistTime + i * refineDt));
            vector2D posAtT = relPos + relVel * t + relAcc * (0.5 * t * t);
            double distAtT = posAtT.norm();

            if (distAtT < minDist) {
                minDist = distAtT;
                minDistTime = t;
            }

            if (distAtT < collisionRadius && t < collisionTime) {
                foundCollision = true;
                collisionTime = t;
            }
        }
    }

    info.willCollide = foundCollision;
    info.collisionTime = foundCollision ? collisionTime : dt;
    info.minDistance = minDist;

    return info;
}

void checkCollisions(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree,
                     double dt) {

#pragma omp parallel
    {
        std::vector<Particle *> query_particles;
        query_particles.reserve(MAX_CAPACITY * 5);
        Bounds query_bounds;

#pragma omp for schedule(dynamic, CHUNK_SIZE)
        for (int i = 0; i < static_cast<int>(particles.size()); i++) {
            auto &particle = particles[i];

            // Velocity-aware search radius: account for distance particle could travel
            double velocityRange = particle->velocity.norm() * dt;
            double range = 2.0 * particle->radius + velocityRange;

            query_bounds.set_bounds(particle->position.x - range, particle->position.y - range,
                                    range * 2, range * 2);
            tree->query(query_bounds, query_particles);

            for (auto *neighbour : query_particles) {
                // Skip if same, already merged, or lower ID (prevents race condition)
                if (particle->id == neighbour->id || neighbour->id <= particle->id ||
                    neighbour->markForDeletion)
                    continue;

                // Use continuous collision detection
                CollisionInfo collision = predictCollision(particle.get(), neighbour, dt);

                if (collision.willCollide) {
                    // Perform perfectly inelastic collision (merge particles)
                    double total_mass = particle->mass + neighbour->mass;

                    particle->velocity = (neighbour->velocity * neighbour->mass +
                                          particle->velocity * particle->mass) /
                                         total_mass;
                    particle->radius = pow(total_mass / particle->mass, 1. / 3.) * particle->radius;
                    particle->mass = total_mass;
                    neighbour->markForDeletion = true;
                    break;
                }
            }
            query_particles.clear();
        }
    }

    /* remove the merged particles from the collision list */
    particles.erase(
        std::remove_if(particles.begin(), particles.end(),
                       [&](const std::shared_ptr<Particle> &p) { return p->markForDeletion; }),
        particles.end());
}
