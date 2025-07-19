#pragma once

#include "global.h"
#include "vector2D.h"

class Bounds
{
public:
    double xmin, ymin, width, height;

    void set_bounds(double _xmin, double _ymin, double _width, double _height)
    {
        this->xmin = _xmin;
        this->ymin = _ymin;
        this->width = _width;
        this->height = _height;
    }

    constexpr bool intersects(Bounds other)
    {
        return !(getLeft() > other.getRight() || getRight() < other.getLeft() || getTop() < other.getBottom() || getBottom() > other.getTop());
    }

    constexpr double getLeft()
    {
        return xmin;
    }

    constexpr double getRight()
    {
        return xmin + width;
    }

    constexpr double getBottom()
    {
        return ymin;
    }

    constexpr double getTop()
    {
        return ymin + height;
    }

    bool contains(vector2D position)
    {
        return ((position.x >= xmin) && (position.x < xmin + width) && (position.y >= ymin) && (position.y < ymin + height));
    }
};

extern Bounds global_bounds;