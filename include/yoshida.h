#pragma once

#include "global.h"
#include "barneshut.h"

const double w0 = -pow(2.0, 1.0 / 3.0) / (2.0 - pow(2.0, 1.0 / 3.0));
const double w1 = 1.0 / (2.0 - pow(2.0, 1.0 / 3.0));

// Coefficients
const double c1 = w1 / 2.0;
const double c2 = (w0 + w1) / 2.0;
const double c3 = c2;
const double c4 = c1;

const double d1 = w1;
const double d2 = w0;
const double d3 = w1;

void drift(std::vector<Particle> &, double);
void kick(std::vector<Particle> &, double);
void yoshidaStep(std::vector<std::shared_ptr<Particle>> &, QuadTree<Particle> *, double);