#include "quadtree.h"
#include "render.h"
#include "interactions.h"
#include <chrono>

Bounds global_bounds;

int main()
{
    srand(5);
    omp_set_num_threads(8);

    Particle *track_particle = nullptr;
    QuadTree<Particle> tree(-250, -250, 500, 500, 1, nullptr);
    global_bounds.set_bounds(-8, -8, 16, 16);
    std::vector<std::shared_ptr<Particle>> particles;

    std::chrono::high_resolution_clock::time_point start, end;
    float fps;

    std::shared_ptr<Particle> particle = std::make_shared<Particle>(0, 0, 0, 0, 0, true);
    particle->mass = 1;
    particle->radius = 0.005;
    particles.push_back(particle);
    for (std::size_t i = 1; i < 100001; i++)
    {
        double dist, angle, x, y, vx, vy, speed;
        dist = ((double)rand() / RAND_MAX) * 4 + 0.25;
        angle = ((double)rand() / RAND_MAX) * 2 * M_PI;
        x = dist * cos(angle);
        y = dist * sin(angle);
        // angle = static_cast<double>(i) / 50 * 2 * M_PI;
        // x = 0.5 * cos(angle);
        // y = 0.5 * sin(angle);
        speed = sqrt(GRAV_G * particles[0]->mass / dist);
        vx = -y / dist * speed;
        vy = x / dist * speed;
        std::shared_ptr<Particle> particle = std::make_shared<Particle>(x, y, vx, vy, i, false);
        particle->mass = 1e-8;
        particle->radius = 1e-8;
        particles.push_back(particle);
    }

    // angle = static_cast<double>(i) / 50 * 2 * M_PI;
    // x = 0.5 * cos(angle);
    // y = 0.5 * sin(angle);
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
        std::shared_ptr<Particle> planet = std::make_shared<Particle>(x, y, vx, vy, particles.size(), true);
        planet->mass = ((double)rand() / RAND_MAX) * 0.001;
        planet->radius = 0.0005;
        particles.push_back(planet);
    }

    std::vector<std::shared_ptr<Particle>> particlesToRemove;
    std::vector<Particle *> query_particles;
    Bounds query_bounds;
    sf::RenderWindow window;
    window.create(sf::VideoMode(GRID_SIZE * CELL_SIZE, GRID_SIZE * CELL_SIZE), "quadtree");
    window.setFramerateLimit(120);
    bool renderTree = false;
    bool pauseSim = true;
    double view_width = 8;
    vector2D view_center(0, 0);
    double time = 0;
    double dt = 0.05;

    for(auto const& particle: particles) {
        tree.insert(particle);
    }

    while (window.isOpen())
    {
        sf::Event event;

        tree.updateParticles(particlesToRemove);

        for (auto &particle : particlesToRemove)
        {
            tree.insert(particle);
        }
        tree.calculateCOM();

        if (!pauseSim)
        {
            updateParticles(particles, &tree, dt);
            // for(auto &particle : particles) {
            //     particle->position.x += -particle->position.y * dt;
            //     particle->position.y += particle->position.x * dt;
            // }
            time += dt;
        }

        if(track_particle != nullptr) {
            view_center = vector2D(track_particle->position);
        } else {
            view_center = vector2D(0, 0);
        }

        while (window.pollEvent(event))
        {
            // if (event.type == sf::Event::MouseMoved)
            // {
            //     sf::Vector2i pos = sf::Mouse::getPosition(window);
            //     vector2D mouseLocation = invtransform(vector2D(static_cast<double>(pos.x), static_cast<double>(pos.y)));

            //     fprintf(stdout, "%d %d   %.3f %.3f\n", event.mouseButton.x, event.mouseButton.y, mouseLocation.x, mouseLocation.y);

            //     query_bounds.set_bounds(mouseLocation.x - 1, mouseLocation.y - 1, 2, 2);

            //     query_particles.clear();
            //     tree.query(query_bounds, query_particles);
            // }
            if (event.type == sf::Event::MouseButtonPressed)
            {
                sf::Vector2i mouse = sf::Mouse::getPosition(window);
                vector2D mouseLocation = invtransform(vector2D(static_cast<double>(mouse.x), static_cast<double>(mouse.y)));

                query_bounds.set_bounds(mouseLocation.x - 0.05 * view_width, mouseLocation.y - 0.05 * view_width, 0.1 * view_width, 0.1 * view_width);

                query_particles.clear();
                tree.query(query_bounds, query_particles);
                track_particle = query_particles[0];
                fprintf(stdout, "%d %d   %.3f %.3f    %d\n", mouse.x, mouse.y, mouseLocation.x, mouseLocation.y, track_particle->id);
                // fprintf(stdout, "%d %.3f (%.3f, %.3f)\n", query_tree->particles.size(), query_tree->totalMass, query_tree->centerOfMass.x, query_tree->centerOfMass.y);
            }
            if (event.type == sf::Event::KeyPressed)
            {
                if (event.key.code == sf::Keyboard::T)
                {
                    renderTree = !renderTree;
                }
                if (event.key.code == sf::Keyboard::C)
                {
                    track_particle = nullptr;
                }
                if (event.key.code == sf::Keyboard::Space)
                {
                    pauseSim = !pauseSim;
                }
            }
            if (event.type == sf::Event::MouseWheelMoved)
            {
                double mouse_delta, new_width;
                mouse_delta = static_cast<double>(event.mouseWheel.delta);
                view_width = std::max(view_width * (1 - mouse_delta / 5), 0.1);
            }
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
        }
        
        global_bounds.set_bounds(-view_width / 2 + view_center.x, -view_width/ 2 + view_center.y, view_width, view_width);

        window.clear(BACKGROUND_COLOR);
        RenderParticles(particles, &window);
        // RenderTestParticles(query_particles, &window);
        RenderTime(time, &window);
        // RenderBounds(query_bounds, &window);
        if (renderTree)
            RenderTree(&tree, &window);
        window.display();
        particlesToRemove.clear();
        // tree.clear();
    }

    return 0;
}
