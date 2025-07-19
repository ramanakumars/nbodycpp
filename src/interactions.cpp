#include "interactions.h"

transport_type TRANSPORT_TYPE = transport_type::YOSHIDA;

void updateParticles(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree, double dt)
{
    transportStep(particles, tree, dt);

    checkCollisions(particles, tree);

    double total_mass = 0, com_x = 0, com_y = 0;
#pragma omp parallel for reduction(+ : com_x, com_y, total_mass) schedule(static, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++)
    {
        if (tree->bounds.contains(particles[i]->position))
        {
            com_x += particles[i]->position.x * particles[i]->mass;
            com_y += particles[i]->position.y * particles[i]->mass;
            total_mass += particles[i]->mass;
        }
    }

    com_x /= total_mass;
    com_y /= total_mass;

    vector2D new_com = {com_x, com_y};
#pragma omp parallel for schedule(static, CHUNK_SIZE)
    for (auto &particle : particles)
    {
        particle->position -= new_com;
    }
}

void transportStep(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree, double dt)
{
    if(TRANSPORT_TYPE == YOSHIDA) {
        yoshidaStep(particles, tree, dt);
    } else if(TRANSPORT_TYPE == LEAPFROG) {
        leapFrogStep(particles, tree, dt);
    } else {
        fprintf(stderr, "%d not a valid transport type\n", TRANSPORT_TYPE);
        exit(1);
    }
}

void checkCollisions(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree)
{

#pragma omp parallel
    {
        std::vector<Particle *> query_particles;
        query_particles.reserve(MAX_CAPACITY * 5);
        Bounds query_bounds;

/* find the colliding particles */
// for (auto &particle : particles)
// {
#pragma omp for schedule(dynamic, CHUNK_SIZE)
        for (int i = 0; i < static_cast<int>(particles.size()); i++)
        {
            auto &particle = particles[i];
            double range = 2.0 * particle->radius;
            query_bounds.set_bounds(particle->position.x - range, particle->position.y - range, range * 2, range * 2);
            tree->query(query_bounds, query_particles);

            for (auto *neighbour : query_particles)
            {
                // Skip if same or already merged
                if (particle->id == neighbour->id || neighbour->markForDeletion)
                    continue;

                if (particle->position.distance(neighbour->position) < particle->radius)
                {
                    double total_mass = particle->mass + neighbour->mass;

                    particle->velocity = (neighbour->velocity * neighbour->mass + particle->velocity * particle->mass) / total_mass;
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
        std::remove_if(
            particles.begin(), particles.end(),
            [&](const std::shared_ptr<Particle> &p)
            {
                return p->markForDeletion;
            }),
        particles.end());
}
