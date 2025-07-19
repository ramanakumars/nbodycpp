#include "barneshut.h"

void getAcceleration(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree)
{
#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++)
    {
        particles[i]->acceleration.zero();
        BarnesHutForce(particles[i].get(), tree, 0.05);
    }
}

vector2D force(const Particle *p1, const Particle *p2)
{
    vector2D diff = p1->position - p2->position;
    double dist = std::max(diff.norm(), 2 * p1->radius); // Prevents singularity
    double invDistCubed = 1.0 / (dist * dist * dist);
    double scale = -GRAV_G * p2->mass * invDistCubed;

    return diff * scale;
}

void BarnesHutForce(Particle *p, const QuadTree<Particle> *tree, double theta)
{
    vector2D diff = p->position - tree->centerOfMass;
    double dist = std::max(diff.norm(), 2 * p->radius);
    double s = tree->bounds.width; // or max(width, height)
    double dist_eff = dist * theta * tree->thetaScale;

    if (s < dist_eff)
    {
        // Acceptable approximation — treat the whole cell as a distant mass
        double inv_r3 = 1.0 / (dist * dist * dist);
        p->acceleration += diff * (-GRAV_G * tree->totalMass * inv_r3);
    }
    else
    {
        if (tree->is_divided)
        {
            // Too close — recurse into children
            for (auto &child : tree->children)
            {
                BarnesHutForce(p, child, theta);
            }
        }
        else
        {
            for (auto const &particle : tree->particles)
            {
                if (p->id != particle->id)
                    p->acceleration += force(p, particle.get());
            }
        }
    }
}
