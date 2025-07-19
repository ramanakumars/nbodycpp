#pragma once

#include "global.h"

class vector2D
{
public:
    double x, y;
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

    void zero()
    {
        x = 0;
        y = 0;
    }

    double norm() const
    {
        return std::sqrt(x * x + y * y);
    }

    double distance(const vector2D b)
    {
        double xx = this->x - b.x;
        double yy = this->y - b.y;
        return std::sqrt(xx * xx + yy * yy);
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