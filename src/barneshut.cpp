/**
 * @file barneshut.cpp
 * @brief Implementation of Barnes-Hut hierarchical force calculation
 */

#include "barneshut.h"

void getAcceleration(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree) {
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++) {
        particles[i]->acceleration.zero();
        BarnesHutForceAndJerk(particles[i].get(), tree, 0.05);
    }
}

/**
 * @brief Compute gravitational acceleration and jerk for particle pair
 *
 * @details Direct particle-particle calculation with softening.
 * Computes both the gravitational force and its time derivative.
 *
 * The jerk formula comes from differentiating Newton's law:
 * d/dt[-Gm*r/r³] = -Gm*[v/r³ - 3(r·v)r/r⁵]
 *
 * @param p1 Particle being acted upon
 * @param p2 Source particle
 * @param[out] acc_out Acceleration output
 * @param[out] jerk_out Jerk output
 */
void forceAndJerk(const Particle *p1, const Particle *p2, vector2D &acc_out, vector2D &jerk_out) {
    vector2D rij = p1->position - p2->position;
    vector2D vij = p1->velocity - p2->velocity;

    double r2 = rij.dot(rij);
    double r = std::sqrt(r2);

    // Softening to prevent singularities
    double r_soft = std::max(r, p1->radius + p2->radius);
    double r_soft2 = r_soft * r_soft;
    double r_soft3 = r_soft2 * r_soft;
    double r_soft5 = r_soft3 * r_soft2;

    // Acceleration: a = -G*m*r/|r|³
    double acc_mag = -GRAV_G * p2->mass / r_soft3;
    acc_out = rij * acc_mag;

    // Jerk: da/dt = -G*m * [v/r³ - 3*(r·v)*r/r⁵]
    // This is the time derivative of acceleration
    double rv_dot = rij.dot(vij);
    double jerk_mag = -GRAV_G * p2->mass / r_soft3;

    jerk_out = vij * jerk_mag - rij * (3.0 * jerk_mag * rv_dot / r_soft2);
}

/**
 * @brief Barnes-Hut tree walk with jerk calculation
 *
 * @details Recursively traverses the quadtree, using the multipole approximation
 * when cells are sufficiently distant. Calculates both acceleration and jerk.
 *
 * Three cases:
 * 1. Cell far enough (s < d*θ): Use multipole approximation
 * 2. Cell subdivided: Recurse into children
 * 3. Leaf node: Direct particle-particle calculation
 *
 * @param p Target particle
 * @param tree Current tree node
 * @param theta Opening angle (accuracy parameter)
 */
void BarnesHutForceAndJerk(Particle *p, const QuadTree<Particle> *tree, double theta) {
    vector2D diff = p->position - tree->centerOfMass;
    double dist = std::max(diff.norm(), 2 * p->radius);
    double s = tree->bounds.width; // or max(width, height)
    double dist_eff = dist * theta * tree->thetaScale;

    if (s < dist_eff) {
        // Acceptable approximation — treat the whole cell as a distant mass
        double r2 = dist * dist;
        double r3 = r2 * dist;

        // Acceleration
        double acc_mag = -GRAV_G * tree->totalMass / r3;
        p->acceleration += diff * acc_mag;

        // Jerk (assuming center of mass velocity ≈ 0 for far-field approximation)
        // For better accuracy, would need to track COM velocity in tree
        // Simplified: jerk ≈ -G*M * 3*(r·v)*r/r⁵
        double rv_dot = diff.dot(p->velocity);
        p->jerk -= diff * (3.0 * acc_mag * rv_dot / r2);
    } else {
        if (tree->is_divided) {
            // Too close — recurse into children
            for (auto &child : tree->children) {
                BarnesHutForceAndJerk(p, child, theta);
            }
        } else {
            for (auto const &particle : tree->particles) {
                if (p->id != particle->id) {
                    // p->acceleration += force(p, particle.get());
                    vector2D acc_contrib, jerk_contrib;
                    forceAndJerk(p, particle.get(), acc_contrib, jerk_contrib);
                    p->acceleration += acc_contrib;
                    p->jerk += jerk_contrib;
                }
            }
        }
    }
}
