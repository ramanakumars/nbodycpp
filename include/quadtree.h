/**
 * @file quadtree.h
 * @brief QuadTree spatial data structure for Barnes-Hut algorithm
 *
 * Implements a hierarchical spatial partitioning structure that enables
 * O(N log N) force calculation in N-body simulations.
 */

#pragma once

#include "global.h"
#include "bounds.h"

/// @brief Maximum particles per leaf node before subdivision
#define MAX_CAPACITY 50

/// @brief Maximum tree depth to prevent infinite recursion
#define MAX_DEPTH 15

/**
 * @class QuadTree
 * @brief Hierarchical spatial partitioning tree for 2D N-body simulation
 *
 * @tparam T Particle type (must have position, mass, velocity members)
 *
 * @details The QuadTree recursively subdivides 2D space into quadrants,
 * storing particles in leaf nodes. This enables:
 * - Fast spatial queries (collision detection, neighbor search)
 * - Barnes-Hut multipole approximation for gravity
 * - O(N log N) force calculation instead of O(N²)
 *
 * Tree structure:
 * - Each node covers a rectangular region (bounds)
 * - Leaf nodes store up to MAX_CAPACITY particles
 * - Internal nodes have 4 children (quadrants: NW, NE, SW, SE)
 * - Subdivision stops at MAX_DEPTH
 *
 * Barnes-Hut properties:
 * - Each node stores total mass and center of mass
 * - Distant groups of particles treated as single mass
 * - Opening angle criterion: s/d < θ
 */
template <class T> class QuadTree {
  public:
    Bounds bounds;                         ///< Spatial region covered by this node
    double totalMass;                      ///< Total mass of all particles in subtree
    double thetaScale;                     ///< Theta scaling factor: (M_ref/M)^alpha
    vector2D centerOfMass;                 ///< Center of mass of all particles in subtree
    int depth;                             ///< Depth in tree (root = 1)
    bool is_divided = false;               ///< True if node is subdivided into children
    std::vector<std::shared_ptr<T>> particles; ///< Particles in this leaf (empty if divided)
    std::array<QuadTree<T> *, 4> children; ///< Child nodes [NW, NE, SW, SE]
    QuadTree<T> *parent;                   ///< Parent node (nullptr for root)

    /**
     * @brief Construct a QuadTree node
     *
     * @param xmin Left edge of region
     * @param ymin Bottom edge of region
     * @param width Width of region
     * @param height Height of region
     * @param _depth Depth in tree (root = 1)
     * @param _parent Pointer to parent node (nullptr for root)
     */
    QuadTree(double xmin, double ymin, double width, double height, int _depth,
                QuadTree *_parent) {
        bounds.set_bounds(xmin, ymin, width, height);
        depth = _depth;
        particles.reserve(MAX_CAPACITY);
        totalMass = 0;
        thetaScale = 0;
        parent = _parent;
    }

    /**
     * @brief Insert a particle into the tree
     *
     * @details Recursively finds the appropriate leaf node for the particle:
     * 1. Check if particle is within bounds
     * 2. If leaf has capacity, add particle
     * 3. If leaf is full, subdivide and redistribute
     * 4. If internal node, recurse to appropriate child
     *
     * @param particle Particle to insert
     * @return True if insertion successful, false if out of bounds
     */
    bool insert(std::shared_ptr<T> particle) {
        if (!bounds.contains(particle->position))
            return false;

        if (((particles.size() < MAX_CAPACITY) && (!is_divided)) || (depth == MAX_DEPTH)) {
            particles.emplace_back(particle);
            return true;
        }
        if (!is_divided) {
            subdivide();
        }
        for (auto &child : children) {
            if (child->insert(particle)) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief Find leaf node intersecting query region
     *
     * @param query_bounds Region to search
     * @param[out] query_tree Pointer to intersecting leaf node
     *
     * @note Returns first intersecting leaf found
     */
    void query_tree(Bounds query_bounds, QuadTree<T> *&query_tree) {
        if (!this->bounds.intersects(query_bounds))
            return;
        if (is_divided) {
            for (auto const &child : children) {
                child->query_tree(query_bounds, query_tree);
            }
        } else {
            query_tree = this;
            return;
        }
    }

    /**
     * @brief Find all particles within query region
     *
     * @details Recursively searches tree, collecting particles that
     * fall within the query bounds.
     *
     * @param query_bounds Region to search
     * @param[out] query_particles Vector to append results to
     *
     * @note Does not clear query_particles - appends results
     * @note Typical use for collision detection, neighbor search
     */
    void query(Bounds query_bounds, std::vector<T *> &query_particles) {
        if (!this->bounds.intersects(query_bounds))
            return;
        if (is_divided) {
            for (auto const &child : children) {
                child->query(query_bounds, query_particles);
            }
        } else {
            for (auto &particle : this->particles) {
                if (query_bounds.contains(particle->position)) {
                    query_particles.emplace_back(particle.get());
                }
            }
        }
    }

    /**
     * @brief Merge child nodes back into parent
     *
     * @details Combines all particles from children into this node
     * and deletes children. Only succeeds if all children are leaf nodes.
     *
     * @return True if merge successful, false otherwise
     *
     * @note Used to coarsen tree when particles leave regions
     */
    bool merge() {
        if (!is_divided)
            return false;

        // Only merge if all children are leaf nodes
        for (auto &child : children) {
            if (child->is_divided)
                return false;
        }

        for (auto &child : children) {
            std::move(child->particles.begin(), child->particles.end(),
                      std::back_inserter(this->particles));
            delete child;
        }
        this->is_divided = false;
        return true;
    }

    /**
     * @brief Calculate center of mass and total mass for subtree
     *
     * @details Recursively computes mass properties for Barnes-Hut algorithm:
     * - For internal nodes: weighted average of children's COM
     * - For leaf nodes: weighted average of particles
     * - Also computes theta scale factor: (M_ref/M)^alpha
     *
     * @note Must be called after particle positions change
     * @note Required before force calculation
     */
    void calculateCOM() {
        centerOfMass = {0, 0};
        totalMass = 0;

        if (is_divided) {
            for (const auto &child : children) {
                child->calculateCOM();
                totalMass += child->totalMass;
                centerOfMass += (child->centerOfMass * child->totalMass);
            }
            centerOfMass /= totalMass;
        } else {
            for (auto &particle : particles) {
                centerOfMass = (centerOfMass * totalMass + particle->position * particle->mass) /
                               (totalMass + particle->mass);
                totalMass += particle->mass;
            }
        }

        thetaScale = std::pow(MASS_REF / totalMass, ALPHA);
    }

    /**
     * @brief Conditionally merge children if particle count is low
     *
     * @details Checks if total particles in all children is below MAX_CAPACITY.
     * If so, and all children are leaves, merges them into this node.
     *
     * @note Automatic tree coarsening for efficiency
     */
    void mergeIfNeeded() {
        if (!is_divided)
            return;
        size_t total_particles = 0;

        for (auto &child : children) {
            total_particles += child->particles.size();
            if (child->is_divided)
                return; // At least one child still subdivided → don't merge
        }

        if (total_particles < MAX_CAPACITY)
            merge();
    }

    /**
     * @brief Remove particles that have left their cells
     *
     * @details Recursively checks if particles are still within bounds.
     * Particles outside bounds are removed and added to particlesToRemove.
     * Automatically merges nodes if they become underpopulated.
     *
     * @param[out] particlesToRemove Vector to collect displaced particles
     *
     * @note Caller must reinsert displaced particles
     * @note Maintains tree consistency during particle motion
     */
    void updateParticles(std::vector<std::shared_ptr<T>> &particlesToRemove) {
        if (is_divided) {
            for (auto const &child : children) {
                child->updateParticles(particlesToRemove);
            }
            mergeIfNeeded();
        } else {
            for (auto it = particles.begin(); it != particles.end();) {
                if (!bounds.contains((*it)->position)) {
                    particlesToRemove.emplace_back(*it);
                    it = particles.erase(it);
                } else {
                    ++it;
                }
            }
        }
    }

  private:
    /**
     * @brief Subdivide node into 4 children
     *
     * @details Creates 4 child quadrants (NW, NE, SW, SE) and redistributes
     * particles from this node to children.
     *
     * Layout:
     * ```
     * [0:NW] [1:NE]
     * [2:SW] [3:SE]
     * ```
     *
     * @note Particles that don't fit in any child remain in parent
     */
    void subdivide() {
        children[0] = new QuadTree<T>(bounds.xmin, bounds.ymin + bounds.height / 2,
                                      bounds.width / 2, bounds.height / 2, depth + 1, this);
        children[1] =
            new QuadTree<T>(bounds.xmin + bounds.width / 2, bounds.ymin + bounds.height / 2,
                            bounds.width / 2, bounds.height / 2, depth + 1, this);
        children[2] = new QuadTree<T>(bounds.xmin, bounds.ymin, bounds.width / 2, bounds.height / 2,
                                      depth + 1, this);
        children[3] = new QuadTree<T>(bounds.xmin + bounds.width / 2, bounds.ymin, bounds.width / 2,
                                      bounds.height / 2, depth + 1, this);
        is_divided = true;

        for (auto it = particles.begin(); it != particles.end();) {
            bool inserted = false;
            for (auto &child : children) {
                if (child->insert(*it)) {
                    inserted = true;
                    break;
                }
            }
            if (inserted) {
                it = particles.erase(it);
            } else {
                ++it;
            }
        }
    }
};
