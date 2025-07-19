#pragma once

#include "global.h"
#include "quadtree.h"
#include "particle.h"

void getAcceleration(std::vector<std::shared_ptr<Particle>> &, QuadTree<Particle> *);

vector2D force(const Particle *, const Particle *);
void BarnesHutForce(Particle *, const QuadTree<Particle> *, double);