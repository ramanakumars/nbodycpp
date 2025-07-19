#include "leapfrog.h"

void leapFrogStep(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree, double dt) {
    getAcceleration(particles, tree);

#pragma omp parallel for schedule(dynamic, CHUNK_SIZE)
    for (int i = 0; i < static_cast<int>(particles.size()); i++)
    {
        auto &p = particles[i];
        Particle temp = *p;

        // First half RK2 step
        temp.position += p->velocity * dt + p->acceleration * (0.5 * dt * dt);

        // Get intermediate acceleration
        temp.acceleration.zero();
        BarnesHutForce(&temp, tree, 0.05);

        // Second half RK2 step
        p->position = temp.position;
        p->velocity += (temp.acceleration + p->acceleration) * (0.5 * dt);
        p->acceleration = temp.acceleration;
    }
}