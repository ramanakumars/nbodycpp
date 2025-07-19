#include "quadtree.h"
#include "render.h"

#define PRIMARY_PARTICLE true
#define NOT_PRIMARY_PARTICLE false
Bounds global_bounds;

int main()
{
    srand(5);
    omp_set_num_threads(8);

    QuadTree<Particle> tree(-250, -250, 500, 500, 1, nullptr);
    global_bounds.set_bounds(-8, -8, 16, 16);
    std::vector<std::shared_ptr<Particle>> particles;

    std::shared_ptr<Particle> particle = std::make_shared<Particle>(0, 0, 0, 0, 0, PRIMARY_PARTICLE);
    particle->mass = 1;
    particle->radius = 0.005;
    particles.push_back(particle);
    
    for (int i = 0; i < 5; i++)
    {
        double dist, angle, x, y, vx, vy, speed;
        dist = ((double)rand() / RAND_MAX) * 5.5 + 0.5;
        angle = ((double)rand() / RAND_MAX) * 2 * M_PI;
        x = dist * cos(angle);
        y = dist * sin(angle);
        speed = sqrt(GRAV_G * particles[0]->mass / dist);
        vx = -y / dist * speed;
        vy = x / dist * speed;
        std::shared_ptr<Particle> planet = std::make_shared<Particle>(x, y, vx, vy, particles.size(), PRIMARY_PARTICLE);
        planet->mass = ((double)rand() / RAND_MAX) * 0.001;
        planet->radius = 0.0005;
        particles.push_back(planet);
    }

    for (std::size_t i = 1; i < 100001; i++)
    {
        double dist, angle, x, y, vx, vy, speed;
        dist = ((double)rand() / RAND_MAX) * 4 + 0.25;
        angle = ((double)rand() / RAND_MAX) * 2 * M_PI;
        x = dist * cos(angle);
        y = dist * sin(angle);
        speed = sqrt(GRAV_G * particles[0]->mass / dist);
        vx = -y / dist * speed;
        vy = x / dist * speed;
        std::shared_ptr<Particle> particle = std::make_shared<Particle>(x, y, vx, vy, particles.size(), NOT_PRIMARY_PARTICLE);
        particle->mass = 1e-8;
        particle->radius = 1e-8;
        particles.push_back(particle);
    }

    for (auto const &particle : particles)
    {
        tree.insert(particle);
    }

    Render renderer = Render(particles, &tree);
    renderer.run(0.01);

    return 0;
}
