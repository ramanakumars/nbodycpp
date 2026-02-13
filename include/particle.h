/**
 * @file particle.h
 * @brief Particle class for N-body simulation
 */

#pragma once

#include "global.h"
#include "vector2D.h"
#include "bounds.h"

/**
 * @class Particle
 * @brief Represents a particle in the N-body simulation
 *
 * @details Stores kinematic state (position, velocity, acceleration, jerk)
 * and physical properties (mass, radius) for gravitational N-body integration.
 */
class Particle
{
public:
    vector2D position;      ///< Current position (x, y)
    vector2D velocity;      ///< Current velocity (vx, vy)
    vector2D acceleration;  ///< Current acceleration (ax, ay)
    vector2D jerk;          ///< Time derivative of acceleration (for Hermite integrator)

    /// @name Predictor variables
    /// @{
    vector2D position_pred; ///< Predicted position (Hermite predictor step)
    vector2D velocity_pred; ///< Predicted velocity (Hermite predictor step)
    /// @}

    int id;                 ///< Unique particle identifier
    double mass;            ///< Particle mass
    double radius;          ///< Particle radius (for collisions and softening)
    bool isPrimary = false; ///< True if particle should be rendered specially
    bool markForDeletion = false; ///< True if particle should be removed (merged)

    /**
     * @brief Construct a new Particle
     *
     * @param x Initial x position
     * @param y Initial y position
     * @param vx Initial x velocity
     * @param vy Initial y velocity
     * @param id Unique particle identifier
     * @param isPrimary Whether this is a primary particle (for rendering)
     *
     * @note Acceleration and jerk are initialized to zero
     */
    Particle(double x, double y, double vx, double vy, int id, bool isPrimary)
    {
        this->position.x = x;
        this->position.y = y;
        this->id = id;
        this->velocity.x = vx;
        this->velocity.y = vy;
        this->acceleration.x = 0;
        this->acceleration.y = 0;
        this->jerk.x = 0;
        this->jerk.y = 0;
        this->isPrimary = isPrimary;
    }

    /**
     * @brief Default copy constructor
     */
    Particle(const Particle &other) = default;
};