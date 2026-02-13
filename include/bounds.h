/**
 * @file bounds.h
 * @brief Axis-aligned bounding box for spatial queries
 */

#pragma once

#include "global.h"
#include "vector2D.h"

/**
 * @class Bounds
 * @brief Axis-aligned rectangular bounding box
 *
 * @details Used for QuadTree cells and spatial query regions.
 * Provides efficient containment and intersection tests.
 */
class Bounds
{
public:
    double xmin;   ///< Minimum x coordinate (left edge)
    double ymin;   ///< Minimum y coordinate (bottom edge)
    double xmax;   ///< Maximum x coordinate (right edge)
    double ymax;   ///< Maximum y coordinate (top edge)
    double width;  ///< Width of bounding box
    double height; ///< Height of bounding box

    /**
     * @brief Set bounding box dimensions
     *
     * @param _xmin Left edge x-coordinate
     * @param _ymin Bottom edge y-coordinate
     * @param _width Width of box
     * @param _height Height of box
     */
    void set_bounds(double _xmin, double _ymin, double _width, double _height)
    {
        xmin = _xmin;
        ymin = _ymin;
        xmax = _xmin + _width;
        ymax = _ymin + _height;
        width = _width;
        height = _height;
    }

    /// @brief Get left edge x-coordinate
    constexpr double getLeft() const { return xmin; }

    /// @brief Get right edge x-coordinate
    constexpr double getRight() const { return xmax; }

    /// @brief Get bottom edge y-coordinate
    constexpr double getBottom() const { return ymin; }

    /// @brief Get top edge y-coordinate
    constexpr double getTop() const { return ymax; }

    /**
     * @brief Check if point is inside bounding box
     *
     * @param position Point to test
     * @return True if position is within bounds (inclusive left/bottom, exclusive right/top)
     */
    inline constexpr bool contains(const vector2D &position) const
    {
        return ((position.x >= getLeft()) && (position.x < getRight()) &&
                (position.y >= getBottom()) && (position.y < getTop()));
    }

    /**
     * @brief Check if this bounding box intersects another
     *
     * @param other Other bounding box
     * @return True if boxes overlap
     */
    inline constexpr bool intersects(const Bounds &other) const
    {
        return !(getLeft() > other.getRight() || getRight() < other.getLeft() ||
                 getTop() < other.getBottom() || getBottom() > other.getTop());
    }
};

/// @brief Global viewing bounds for rendering
extern Bounds global_bounds;