#pragma once

#include "global.h"
#include "quadtree.h"
#include "particle.h"
#include "yoshida.h"
#include "leapfrog.h"

void transportStep(std::vector<std::shared_ptr<Particle>> &particles, QuadTree<Particle> *tree, double dt);
void updateParticles(std::vector<std::shared_ptr<Particle>> &, QuadTree<Particle> *, double);
void checkCollisions(std::vector<std::shared_ptr<Particle>> &, QuadTree<Particle> *);
