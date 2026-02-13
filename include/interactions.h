/**
 * @file interactions.h
 * @brief Particle interactions including integration and collision handling
 */

#pragma once

#include "global.h"
#include "quadtree.h"
#include "particle.h"

/**
 * @brief Perform one integration timestep using selected integrator
 *
 * @details Dispatches to the appropriate integrator based on TRANSPORT_TYPE:
 * - YOSHIDA: 4th order symplectic integrator
 * - RK2: 2nd order Runge-Kutta
 * - HERMITE: 4th order Hermite predictor-corrector
 *
 * @param particles Particle vector
 * @param tree QuadTree for force calculation
 * @param dt Timestep size
 */
void transportStep(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree, double dt);

/**
 * @brief Main update function: integrate and handle collisions
 *
 * @details Performs one complete timestep:
 * 1. Integrate particle positions/velocities
 * 2. Check and resolve collisions
 * 3. Recenter system to center of mass
 *
 * @param particles Particle vector
 * @param tree QuadTree for force calculation
 * @param dt Timestep size
 */
void updateParticles(std::vector<std::shared_ptr<Particle>> &, QuadTree<Particle> *, double);

/**
 * @brief Detect and resolve particle collisions
 *
 * @details Uses quadtree spatial partitioning and continuous collision detection
 * to find colliding particles. Merges particles via perfectly inelastic collisions,
 * conserving momentum and mass.
 *
 * Features:
 * - Velocity-aware search radius
 * - Continuous collision detection (accounts for gravitational trajectories)
 * - Thread-safe with directional merging (ID-based)
 * - Removes merged particles from simulation
 *
 * @param particles Particle vector (modified in-place)
 * @param tree QuadTree for spatial queries
 * @param dt Timestep (for continuous collision detection)
 *
 * @note Parallelized with OpenMP
 * @note Uses race-condition prevention via ID ordering
 */
void checkCollisions(std::vector<std::shared_ptr<Particle>> &, QuadTree<Particle> *, double dt);

/**
 * @struct CollisionInfo
 * @brief Result of collision prediction between two particles
 */
struct CollisionInfo {
    bool willCollide;       ///< True if particles will collide during timestep
    double collisionTime;   ///< Time within [0, dt] when collision occurs
    double minDistance;     ///< Minimum separation distance during timestep
};

/**
 * @brief Predict if and when two particles will collide
 *
 * @details Uses continuous collision detection with gravitational trajectory
 * approximation. Samples the quadratic trajectory (constant acceleration)
 * to find minimum approach distance and collision time.
 *
 * Trajectory model: r(t) = r₀ + v₀*t + ½*a*t²
 *
 * @param p1 First particle
 * @param p2 Second particle
 * @param dt Timestep to check within
 * @return CollisionInfo structure with prediction results
 *
 * @note Uses constant acceleration approximation
 * @note Includes gravitational effects in trajectory
 * @note Samples trajectory at multiple points for accuracy
 */
CollisionInfo predictCollision(const Particle *p1, const Particle *p2, double dt);
