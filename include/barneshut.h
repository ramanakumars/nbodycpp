/**
 * @file barneshut.h
 * @brief Barnes-Hut algorithm for hierarchical N-body force calculation
 *
 * Implements the Barnes-Hut tree algorithm for computing gravitational forces
 * in O(N log N) time instead of O(N²) for direct summation.
 */

#pragma once

#include "global.h"
#include "quadtree.h"
#include "particle.h"

/**
 * @brief Calculate acceleration for all particles using Barnes-Hut algorithm
 *
 * @param particles Vector of particles to update
 * @param tree QuadTree structure for hierarchical force calculation
 *
 * @note Zeros acceleration before calculation
 * @note Parallelized with OpenMP
 */
void getAcceleration(std::vector<std::shared_ptr<Particle>> &, QuadTree<Particle> *);

/**
 * @brief Calculate gravitational force and jerk between two particles
 *
 * @details Computes both acceleration and its time derivative (jerk) for
 * gravitational interaction between a particle pair.
 *
 * Formulas:
 * - Acceleration: a = -G*m*r/|r|³
 * - Jerk: da/dt = -G*m * [v/r³ - 3*(r·v)*r/r⁵]
 *
 * @param p1 First particle (particle being acted upon)
 * @param p2 Second particle (source of force)
 * @param[out] acc_out Acceleration contribution
 * @param[out] jerk_out Jerk contribution (time derivative of acceleration)
 *
 * @note Uses gravitational softening based on particle radii
 */
void forceAndJerk(const Particle *, const Particle *, vector2D &, vector2D &);

/**
 * @brief Recursively calculate force and jerk using Barnes-Hut tree
 *
 * @details Implements the Barnes-Hut multipole approximation with jerk calculation.
 * Uses opening angle criterion to determine when to approximate distant particles
 * as a single mass.
 *
 * Opening criterion: s/d < θ
 * - s = cell size
 * - d = distance to cell center of mass
 * - θ = opening angle (smaller = more accurate)
 *
 * @param p Particle to calculate forces for
 * @param tree QuadTree node to evaluate
 * @param theta Opening angle parameter (typical: 0.05-0.5)
 *
 * @note Accumulates into p->acceleration and p->jerk
 * @note For far-field, assumes COM velocity ≈ 0 (simplification)
 */
void BarnesHutForceAndJerk(Particle *p, const QuadTree<Particle> *tree,
                          double theta);
