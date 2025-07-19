#include "render.h"

sf::Color BACKGROUND_COLOR(0, 0, 0);
sf::Color PARTICLE_COLOR(255, 255, 255);
sf::Color PRIMARY_COLOR(255, 0, 255);

sf::Vector2<float> transform(vector2D vec)
{
    double scalex = (GRID_SIZE * CELL_SIZE - 2 * START) / (global_bounds.getRight() - global_bounds.getLeft());
    double scaley = (GRID_SIZE * CELL_SIZE - 2 * START) / (global_bounds.getTop() - global_bounds.getBottom());
    double x, y;

    x = (vec.x - global_bounds.getLeft()) * scalex + START;
    y = (vec.y - global_bounds.getBottom()) * scaley + START;

    return sf::Vector2<float>({static_cast<float>(x), static_cast<float>(y)});
}

vector2D invtransform(vector2D vec)
{
    double scalex = (GRID_SIZE * CELL_SIZE - 2 * START) / (global_bounds.getRight() - global_bounds.getLeft());
    double scaley = (GRID_SIZE * CELL_SIZE - 2 * START) / (global_bounds.getTop() - global_bounds.getBottom());
    double x, y;

    x = (vec.x - START) / scalex + global_bounds.getLeft();
    y = (vec.y - START) / scaley + global_bounds.getBottom();

    return vector2D(x, y);
}

void RenderTestParticles(const std::vector<Particle *> particles, sf::RenderWindow *window)
{
    const float size = 2;
    for (auto const &p : particles)
    {
        sf::CircleShape circle(size);

        circle.setFillColor(sf::Color(0, 255, 0));
        circle.setPosition(transform(p->position) - sf::Vector2<float>({static_cast<float>(size) / 2, static_cast<float>(size) / 2}));
        window->draw(circle);
    }
}

void RenderParticles(const std::vector<std::shared_ptr<Particle>> &particles, sf::RenderWindow *window)
{
    double two_mu = 2 * GRAV_G * particles[0]->mass;
    sf::VertexArray points(sf::Points, particles.size());
    std::vector<sf::CircleShape> circles;
    circles.reserve(15);

    for (int i = 0; i < particles.size(); ++i)
    {
        const auto &p = particles[i];
        if (!global_bounds.contains(p->position))
            continue;
        if (p->isPrimary)
        {
            double size = (log10(p->radius) + 5);
            sf::CircleShape circle(size);

            circle.setFillColor(PRIMARY_COLOR);
            circle.setPosition(transform(p->position) - sf::Vector2<float>({static_cast<float>(size) / 2, static_cast<float>(size) / 2}));

            circles.push_back(std::move(circle));
        }
        points[i].position = transform(p->position);
        double dist = (p->position - particles[0]->position).norm();
        if (p->velocity.norm() < sqrt(two_mu / dist))
        {
            points[i].color = PARTICLE_COLOR;
        }
        else
        {
            points[i].color = sf::Color::Red;
        }
    }

    window->draw(points);
    // Serial draw loop
    for (auto &circle : circles)
        window->draw(circle);
}

void RenderTree(QuadTree<Particle> *tree, sf::RenderWindow *window)
{
    if (tree->is_divided)
    {
        for (auto const &child : tree->children)
        {
            RenderTree(child, window);
        }
    } else {
        if(!tree->bounds.intersects(global_bounds)) return;
        RenderBounds(tree->bounds, window);
    }
}

void RenderBounds(Bounds bounds, sf::RenderWindow *window)
{
    RenderLine(
        vector2D(bounds.getLeft(), bounds.getBottom()),
        vector2D(bounds.getRight(), bounds.getBottom()),
        window);
    RenderLine(
        vector2D(bounds.getRight(), bounds.getBottom()),
        vector2D(bounds.getRight(), bounds.getTop()),
        window);
    RenderLine(
        vector2D(bounds.getRight(), bounds.getTop()),
        vector2D(bounds.getLeft(), bounds.getTop()),
        window);
    RenderLine(
        vector2D(bounds.getLeft(), bounds.getTop()),
        vector2D(bounds.getLeft(), bounds.getBottom()),
        window);
}

void RenderLine(vector2D position1, vector2D position2, sf::RenderWindow *window)
{
    std::vector<sf::Vertex> line;
    line.push_back(sf::Vertex{transform(position1)});
    line.push_back(sf::Vertex{transform(position2)});
    window->draw(line.data(), line.size(), sf::PrimitiveType::Lines);
}

void RenderTime(double time, sf::RenderWindow *window)
{
    static int initialized;
    static sf::Font font;
    char message[200];

    if (!initialized)
    {
        if (!font.loadFromFile("../fonts/arial.ttf"))
        {
            std::cerr << "Failed to load font!" << std::endl;
            return;
        }
        initialized = 1;
    }

    sprintf(message, "Time: %.2f years", time);

    sf::Text text;
    text.setFont(font);
    text.setString(message);
    text.setCharacterSize(18);
    text.setFillColor(sf::Color::White);
    sf::FloatRect textRect = text.getLocalBounds();
    text.setOrigin(textRect.width / 2, textRect.height / 2);
    text.setPosition(GRID_SIZE / 2, START / 2);
    window->draw(text);
}