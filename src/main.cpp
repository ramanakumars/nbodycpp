/**
 * @file main.cpp
 * @brief N-body simulation of planetary system with 100k test particles
 *
 * @details Simulates a central star with planets and a debris disk.
 * Uses Barnes-Hut algorithm with quadtree spatial partitioning for
 * efficient gravity calculation.
 *
 * System setup:
 * - 1 central star (mass = 1.0)
 * - 5 planets in circular orbits
 * - 100,000 massless test particles (debris disk)
 *
 * Integration: Hermite 4th order (default), configurable via TRANSPORT_TYPE
 * Timestep: 0.01 time units
 * Visualization: SFML-based real-time rendering
 */

#include "quadtree.h"
#include "render.h"

#define PRIMARY_PARTICLE true       ///< Particle rendered as colored circle
#define NOT_PRIMARY_PARTICLE false  ///< Particle rendered as point
Bounds global_bounds;               ///< Global viewing bounds

/**
 * @brief Main entry point for N-body simulation
 *
 * @details Initialization sequence:
 * 1. Set random seed for reproducibility
 * 2. Configure OpenMP threads
 * 3. Create QuadTree and initialize bounds
 * 4. Generate particle system:
 *    - Central star at origin
 *    - 5 planets in Keplerian orbits (0.5-6 units radius)
 *    - 100k test particles in disk (0.25-4.25 units radius)
 * 5. Insert particles into QuadTree
 * 6. Launch renderer with dt=0.01
 *
 * @return 0 on successful completion
 */
int main()
{
    srand(5);               // Fixed seed for reproducibility
    omp_set_num_threads(8); // Parallel computation with 8 threads

    // Create QuadTree covering [-250, -250] to [250, 250]
    QuadTree<Particle> tree(-250, -250, 500, 500, 1, nullptr);

    // Set initial viewing bounds [-8, -8] to [8, 8]
    global_bounds.set_bounds(-8, -8, 16, 16);

    std::vector<std::shared_ptr<Particle>> particles;

    // Create central star
    std::shared_ptr<Particle> particle = std::make_shared<Particle>(0, 0, 0, 0, 0, PRIMARY_PARTICLE);
    particle->mass = 1;        // Star mass = 1.0 (unit mass)
    particle->radius = 0.005;  // Star radius
    particles.push_back(particle);

    // Create 5 planets in circular Keplerian orbits
    for (int i = 0; i < 5; i++)
    {
        double dist, angle, x, y, vx, vy, speed;
        dist = ((double)rand() / RAND_MAX) * 5.5 + 0.5;  // Random radius [0.5, 6]
        angle = ((double)rand() / RAND_MAX) * 2 * M_PI;   // Random angle
        x = dist * cos(angle);
        y = dist * sin(angle);

        // Circular orbit velocity: v = sqrt(GM/r)
        speed = sqrt(GRAV_G * particles[0]->mass / dist);
        vx = -y / dist * speed;  // Tangential velocity
        vy = x / dist * speed;

        std::shared_ptr<Particle> planet = std::make_shared<Particle>(x, y, vx, vy, particles.size(), PRIMARY_PARTICLE);
        planet->mass = ((double)rand() / RAND_MAX) * 0.001;  // Planet mass [0, 0.001]
        planet->radius = 0.0005;
        particles.push_back(planet);
    }

    // Create 100,000 test particles (debris disk)
    for (std::size_t i = 1; i < 100001; i++)
    {
        double dist, angle, x, y, vx, vy, speed;
        dist = ((double)rand() / RAND_MAX) * 4 + 0.25;   // Random radius [0.25, 4.25]
        angle = ((double)rand() / RAND_MAX) * 2 * M_PI;   // Random angle
        x = dist * cos(angle);
        y = dist * sin(angle);

        // Circular orbit velocity
        speed = sqrt(GRAV_G * particles[0]->mass / dist);
        vx = -y / dist * speed;
        vy = x / dist * speed;

        std::shared_ptr<Particle> particle = std::make_shared<Particle>(x, y, vx, vy, particles.size(), NOT_PRIMARY_PARTICLE);
        particle->mass = 1e-8;    // Nearly massless test particles
        particle->radius = 1e-8;  // Very small radius
        particles.push_back(particle);
    }

    // Insert all particles into QuadTree
    for (auto const &particle : particles)
    {
        tree.insert(particle);
    }

    // Create renderer and run simulation
    Render renderer = Render(particles, &tree);
    renderer.run(0.01);  // Timestep dt = 0.01

    return 0;
}
