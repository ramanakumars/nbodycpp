/**
 * @file vector2D.h
 * @brief 2D vector class for N-body simulation
 */

#pragma once

#include "global.h"

/**
 * @class vector2D
 * @brief Two-dimensional vector for positions, velocities, and accelerations
 *
 * @details Provides vector arithmetic operations and utility functions
 * optimized for physics simulations.
 */
class vector2D
{
public:
    double x; ///< x component
    double y; ///< y component
    // Constructors
    constexpr vector2D() noexcept : x(0), y(0) {}
    constexpr vector2D(double _x, double _y) noexcept : x(_x), y(_y) {}
    constexpr vector2D(const vector2D &other) noexcept = default;
    constexpr vector2D &operator=(const vector2D &other) noexcept = default;

    // Compound operators
    constexpr vector2D &operator+=(const vector2D &other) noexcept
    {
        x += other.x;
        y += other.y;
        return *this;
    }

    constexpr vector2D &operator-=(const vector2D &other) noexcept
    {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    constexpr vector2D &operator*=(const double val) noexcept
    {
        x *= val;
        y *= val;
        return *this;
    }

    constexpr vector2D &operator/=(const double val) noexcept
    {
        x /= val;
        y /= val;
        return *this;
    }

    /**
     * @brief Set vector to zero
     */
    void zero()
    {
        x = 0;
        y = 0;
    }

    /**
     * @brief Calculate Euclidean norm (magnitude) of vector
     * @return ||v|| = sqrt(x² + y²)
     */
    double norm() const
    {
        return std::sqrt(x * x + y * y);
    }

    /**
     * @brief Calculate distance to another vector
     * @param b Other vector
     * @return ||this - b||
     */
    double distance(const vector2D b) const
    {
        double xx = this->x - b.x;
        double yy = this->y - b.y;
        return std::sqrt(xx * xx + yy * yy);
    }

    /**
     * @brief Calculate dot product with another vector
     * @param other Other vector
     * @return this·other = x*other.x + y*other.y
     */
    double dot(const vector2D &other) const
    {
        return x * other.x + y * other.y;
    }
};

// Addition
constexpr vector2D operator+(vector2D a, const vector2D &b)
{
    return a += b;
}

// Subtraction
constexpr vector2D operator-(vector2D a, const vector2D &b)
{
    return a -= b;
}

// Scalar multiplication
constexpr vector2D operator*(vector2D a, const double val)
{
    return a *= val;
}

// Scalar division
constexpr vector2D operator/(vector2D a, const double val)
{
    return a /= val;
}