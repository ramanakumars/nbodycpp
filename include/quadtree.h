#pragma once

#include "global.h"
#include "bounds.h"

#define MAX_CAPACITY 50
#define MAX_DEPTH 15

template <class T> class QuadTree {
  public:
    Bounds bounds;
    double totalMass, thetaScale;
    vector2D centerOfMass;
    int depth;
    bool is_divided = false;
    std::vector<std::shared_ptr<T>> particles;
    std::array<QuadTree<T> *, 4> children;
    QuadTree<T> *parent;

    QuadTree(double xmin, double ymin, double width, double height, int _depth,
                QuadTree *_parent) {
        bounds.set_bounds(xmin, ymin, width, height);
        depth = _depth;
        particles.reserve(MAX_CAPACITY);
        totalMass = 0;
        thetaScale = 0;
        parent = _parent;
    }

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
