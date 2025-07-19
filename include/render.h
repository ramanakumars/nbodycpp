#pragma once

#include "global.h"
#include "particle.h"
#include "bounds.h"
#include "quadtree.h"
#include <SFML/Graphics.hpp>

#define CELL_SIZE 1
#define GRID_SIZE 800
#define START 50

extern sf::Color BACKGROUND_COLOR; 
extern sf::Color PARTICLE_COLOR; 
extern sf::Color PRIMARY_COLOR; 

sf::Vector2<float> transform(vector2D);
vector2D invtransform(vector2D);
void RenderTestParticles(const std::vector<Particle *> particles, sf::RenderWindow *window);
void RenderParticles(const std::vector<std::shared_ptr<Particle>> &, sf::RenderWindow *);
void RenderTree(QuadTree<Particle> *, sf::RenderWindow *);
void RenderBounds(Bounds, sf::RenderWindow *);
void RenderLine(vector2D, vector2D, sf::RenderWindow *);
void RenderTime(double, sf::RenderWindow *);