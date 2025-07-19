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

    void checkBounds(Bounds bounds)
    {
        return;
        if ((position.x > bounds.getRight()) || (position.x < bounds.getLeft()))
        {
            velocity.x = -velocity.x;
            position.x = fmax(bounds.getRight(), position.x);
            position.x = fmin(bounds.getLeft(), position.x);
        }
        if ((position.y > bounds.getTop()) || (position.y < bounds.getBottom()))
        {
            velocity.y = -velocity.y;
            position.y = fmax(bounds.getTop(), position.y);
            position.y = fmin(bounds.getBottom(), position.y);
        }
    }
};