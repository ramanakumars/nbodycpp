/**
 * @file render.h
 * @brief Graphics rendering for N-body simulation
 */

#pragma once

#include "global.h"
#include "particle.h"
#include "bounds.h"
#include "quadtree.h"
#include "interactions.h"
#include <SFML/Graphics.hpp>

#define CELL_SIZE 1   ///< Pixel size of grid cells
#define GRID_SIZE 800 ///< Window size in pixels
#define START 50      ///< Margin from window edge

/// @name Rendering colors
/// @{
extern sf::Color BACKGROUND_COLOR; ///< Background color (black)
extern sf::Color PARTICLE_COLOR;   ///< Normal particle color (white)
extern sf::Color PRIMARY_COLOR;    ///< Primary particle color (magenta)
/// @}

/**
 * @brief Transform simulation coordinates to screen coordinates
 * @param vec Position in simulation space
 * @return Position in screen space
 */
sf::Vector2<float> transform(vector2D vec);

/**
 * @brief Transform screen coordinates to simulation coordinates
 * @param vec Position in screen space
 * @return Position in simulation space
 */
vector2D invtransform(vector2D vec);

/**
 * @class Render
 * @brief Main rendering class for visualization
 *
 * @details Handles SFML window creation, event processing, and rendering
 * of particles, quadtree structure, and UI elements.
 *
 * Features:
 * - Real-time particle visualization
 * - QuadTree structure overlay
 * - Particle tracking (click to follow)
 * - Zoom and pan controls
 * - Velocity-based coloring (bound vs unbound orbits)
 * - Pause/resume simulation
 */
class Render
{
public:
    /**
     * @brief Construct renderer
     *
     * @param _particles Reference to particle vector
     * @param _tree Pointer to QuadTree structure
     */
    Render(std::vector<std::shared_ptr<Particle>> &_particles, QuadTree<Particle> *_tree)
        : particles(_particles), tree(_tree)
    {
        window.create(sf::VideoMode(GRID_SIZE * CELL_SIZE, GRID_SIZE * CELL_SIZE), "quadtree");
        window.setFramerateLimit(120);
    }

    /// @brief Render debug particles (green circles)
    void renderTestParticles(const std::vector<Particle *> particles);

    /// @brief Render all particles with velocity-based coloring
    void renderParticles(const std::vector<std::shared_ptr<Particle>> &);

    /// @brief Render QuadTree structure overlay
    void renderTree(QuadTree<Particle> *);

    /// @brief Render a bounding box outline
    void renderBounds(Bounds);

    /// @brief Render a line between two points
    void renderLine(vector2D, vector2D);

    /// @brief Render simulation time display
    void renderTime(double);

    /// @brief Render information about tracked particle
    void renderParticleInfo(Particle *);

    /**
     * @brief Main render loop
     *
     * @details Handles:
     * - Event processing (mouse, keyboard)
     * - Simulation stepping
     * - Rendering all elements
     * - QuadTree updates
     *
     * Controls:
     * - Mouse click: Track particle
     * - Space: Pause/resume
     * - T: Toggle tree visualization
     * - C: Clear particle tracking
     * - Mouse wheel: Zoom
     *
     * @param dt Timestep for integration
     */
    void run(double);

private:
    sf::RenderWindow window;                               ///< SFML render window
    std::vector<std::shared_ptr<Particle>> &particles;     ///< Particle array reference
    QuadTree<Particle> *tree;                              ///< QuadTree pointer
    Particle *track_particle = nullptr;                    ///< Currently tracked particle
    std::vector<Particle *> query_particles;               ///< Query result buffer
    Bounds query_bounds;                                   ///< Query region
    bool shouldRenderTree = false;                         ///< Show QuadTree overlay
    bool pauseSim = true;                                  ///< Simulation paused
    double view_width = 8;                                 ///< View width in sim units
    vector2D view_center;                                  ///< View center position

    /**
     * @brief Find particle at mouse cursor position
     *
     * @details Uses QuadTree spatial query to find nearest particle
     * to mouse click location. Updates track_particle.
     */
    void findTrackParticle();
};
