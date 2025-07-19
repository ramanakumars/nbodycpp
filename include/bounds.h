#pragma once

#include "global.h"
#include "vector2D.h"

class Bounds
{
public:
    double xmin, ymin, xmax, ymax, width, height;

    void set_bounds(double _xmin, double _ymin, double _width, double _height)
    {
        xmin = _xmin;
        ymin = _ymin;
        xmax = _xmin + _width;
        ymax = _ymin + _height;
        width = _width;
        height = _height;
    }

    constexpr double getLeft() const
    {
        return xmin;
    }

    constexpr double getRight() const
    {
        return xmax;
    }

    constexpr double getBottom() const
    {
        return ymin;
    }

    constexpr double getTop() const
    {
        return ymax;
    }

    inline constexpr bool contains(const vector2D &position) const
    {
        return ((position.x >= getLeft()) && (position.x < getRight()) && (position.y >= getBottom()) && (position.y < getTop()));
    }

    inline constexpr bool intersects(const Bounds &other) const
    {
        return !(getLeft() > other.getRight() || getRight() < other.getLeft() || getTop() < other.getBottom() || getBottom() > other.getTop());
    }
};

extern Bounds global_bounds;