#pragma once

#include "global.h"
#include "vector2D.h"
#include "bounds.h"

class Particle
{
public:
    vector2D position, velocity, acceleration;
    int id;
    double mass, radius;
    bool isPrimary = false;
    bool markForDeletion = false;

    Particle(double x, double y, double vx, double vy, int id, bool isPrimary)
    {
        this->position.x = x;
        this->position.y = y;
        this->id = id;
        this->velocity.x = vx;
        this->velocity.y = vy;
        this->acceleration.x = 0;
        this->acceleration.y = 0;
        this->isPrimary = isPrimary;
    }

    Particle(const Particle &other) = default;
};