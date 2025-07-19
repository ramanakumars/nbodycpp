#pragma once

#include "global.h"
#include "particle.h"
#include "bounds.h"
#include "quadtree.h"
#include "interactions.h"
#include <SFML/Graphics.hpp>

#define CELL_SIZE 1
#define GRID_SIZE 800
#define START 50

extern sf::Color BACKGROUND_COLOR;
extern sf::Color PARTICLE_COLOR;
extern sf::Color PRIMARY_COLOR;

sf::Vector2<float> transform(vector2D);
vector2D invtransform(vector2D);

class Render
{
public:
    Render(std::vector<std::shared_ptr<Particle>> &_particles, QuadTree<Particle> *_tree) : particles(_particles),
                                                                                            tree(_tree)
    {
        window.create(sf::VideoMode(GRID_SIZE * CELL_SIZE, GRID_SIZE * CELL_SIZE), "quadtree");
        window.setFramerateLimit(120);
    }

    void renderTestParticles(const std::vector<Particle *> particles);
    void renderParticles(const std::vector<std::shared_ptr<Particle>> &);
    void renderTree(QuadTree<Particle> *);
    void renderBounds(Bounds);
    void renderLine(vector2D, vector2D);
    void renderTime(double);
    void renderParticleInfo(Particle *);
    void run(double);

private:
    sf::RenderWindow window;
    std::vector<std::shared_ptr<Particle>> &particles;
    QuadTree<Particle> *tree;
    Particle *track_particle = nullptr;
    std::vector<Particle *> query_particles;
    Bounds query_bounds;
    bool shouldRenderTree = false;
    bool pauseSim = true;
    double view_width = 8;
    vector2D view_center;

    void findTrackParticle();

};
