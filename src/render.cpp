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

void Render::renderTestParticles(const std::vector<Particle *> particles)
{
    const float size = 2;
    for (auto const &p : particles)
    {
        sf::CircleShape circle(size);

        circle.setFillColor(sf::Color(0, 255, 0));
        circle.setPosition(transform(p->position) - sf::Vector2<float>({static_cast<float>(size) / 2, static_cast<float>(size) / 2}));
        window.draw(circle);
    }
}

void Render::renderParticles(const std::vector<std::shared_ptr<Particle>> &particles)
{
    Particle *central_particle;
    track_particle? central_particle = track_particle : central_particle = particles[0].get();
    double two_mu = 2 * GRAV_G * central_particle->mass;
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
        else
        {
            points[i].position = transform(p->position);
            double dist = (p->position - central_particle->position).norm();
            if ((p->velocity - central_particle->velocity).norm() < sqrt(two_mu / dist))
            {
                points[i].color = PARTICLE_COLOR;
            }
            else
            {
                points[i].color = sf::Color::Red;
            }
        }
    }

    window.draw(points);
    // Serial draw loop
    for (auto &circle : circles)
        window.draw(circle);
}

void Render::renderTree(QuadTree<Particle> *tree)
{
    if (tree->is_divided)
    {
        for (auto const &child : tree->children)
        {
            renderTree(child);
        }
    }
    else
    {
        if (!tree->bounds.intersects(global_bounds))
            return;
        renderBounds(tree->bounds);
    }
}

void Render::renderBounds(Bounds bounds)
{
    renderLine(
        vector2D(bounds.getLeft(), bounds.getBottom()),
        vector2D(bounds.getRight(), bounds.getBottom()));
    renderLine(
        vector2D(bounds.getRight(), bounds.getBottom()),
        vector2D(bounds.getRight(), bounds.getTop()));
    renderLine(
        vector2D(bounds.getRight(), bounds.getTop()),
        vector2D(bounds.getLeft(), bounds.getTop()));
    renderLine(
        vector2D(bounds.getLeft(), bounds.getTop()),
        vector2D(bounds.getLeft(), bounds.getBottom()));
}

void Render::renderLine(vector2D position1, vector2D position2)
{
    std::vector<sf::Vertex> line;
    line.push_back(sf::Vertex{transform(position1)});
    line.push_back(sf::Vertex{transform(position2)});
    window.draw(line.data(), line.size(), sf::PrimitiveType::Lines);
}

void Render::renderTime(double time)
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
    window.draw(text);
}

void Render::renderParticleInfo(Particle *particle)
{
    static int initialized;
    static sf::Font font;
    char id_message[200], position_message[200], velocity_message[200];

    if (!initialized)
    {
        if (!font.loadFromFile("../fonts/arial.ttf"))
        {
            std::cerr << "Failed to load font!" << std::endl;
            return;
        }
        initialized = 1;
    }

    sprintf(id_message, "Particle ID: %6d  Mass: %5.3e", particle->id, particle->mass);
    sprintf(position_message, "Position: (%5.2f, %5.2f)", particle->position.x, particle->position.y);
    sprintf(velocity_message, "Velocity: (%5.2f, %5.2f) %5.2f", particle->velocity.x, particle->velocity.y, particle->velocity.norm());

    sf::Text id_text, position_text, velocity_text;
    id_text.setFont(font);
    id_text.setString(id_message);
    id_text.setCharacterSize(14);
    id_text.setFillColor(sf::Color::White);
    sf::FloatRect id_textRect = id_text.getLocalBounds();
    id_text.setOrigin(0, id_textRect.height / 2);
    id_text.setPosition(GRID_SIZE * 2 / 3, START / 2 - 14);

    position_text.setFont(font);
    position_text.setString(position_message);
    position_text.setCharacterSize(14);
    position_text.setFillColor(sf::Color::White);
    sf::FloatRect position_textRect = position_text.getLocalBounds();
    position_text.setOrigin(0, position_textRect.height / 2);
    position_text.setPosition(GRID_SIZE * 2 / 3, START / 2);

    velocity_text.setFont(font);
    velocity_text.setString(velocity_message);
    velocity_text.setCharacterSize(14);
    velocity_text.setFillColor(sf::Color::White);
    sf::FloatRect velocity_textRect = velocity_text.getLocalBounds();
    velocity_text.setOrigin(0, velocity_textRect.height / 2);
    velocity_text.setPosition(GRID_SIZE * 2 / 3, START / 2 + 14);

    window.draw(id_text);
    window.draw(position_text);
    window.draw(velocity_text);
}

void Render::run(double dt)
{
    double time = 0;
    view_center = {0, 0};

    while (window.isOpen())
    {
        sf::Event event;

        std::vector<std::shared_ptr<Particle>> particlesToRemove;
        particlesToRemove.reserve(10000);

        tree->updateParticles(particlesToRemove);

        for (auto &particle : particlesToRemove)
        {
            tree->insert(particle);
        }
        tree->calculateCOM();

        if (!pauseSim)
        {
            updateParticles(particles, tree, dt);
            time += dt;
        }

        if (track_particle != nullptr)
        {
            view_center = track_particle->position;
        }
        else
        {
            view_center = {0, 0};
        }

        while (window.pollEvent(event))
        {
            // if (event.type == sf::Event::MouseMoved)
            // {
            //     sf::Vector2i pos = sf::Mouse::getPosition(window);
            //     vector2D mouseLocation = invtransform(vector2D(static_cast<double>(pos.x), static_cast<double>(pos.y)));

            //     fprintf(stdout, "%d %d   %.3f %.3f\n", event.mouseButton.x, event.mouseButton.y, mouseLocation.x, mouseLocation.y);

            //     query_bounds.set_bounds(mouseLocation.x - 1, mouseLocation.y - 1, 2, 2);

            //     query_particles.clear();
            //     tree.query(query_bounds, query_particles);
            // }
            if (event.type == sf::Event::MouseButtonPressed)
            {
                findTrackParticle();
            }

            if (event.type == sf::Event::KeyPressed)
            {
                if (event.key.code == sf::Keyboard::T)
                    shouldRenderTree = !shouldRenderTree;
                if (event.key.code == sf::Keyboard::C)
                    track_particle = nullptr;
                if (event.key.code == sf::Keyboard::Space)
                    pauseSim = !pauseSim;
            }

            if (event.type == sf::Event::MouseWheelMoved)
            {
                double mouse_delta;
                mouse_delta = static_cast<double>(event.mouseWheel.delta);
                view_width = std::max(view_width * (1 - mouse_delta / 5), 0.1);
            }
            if (event.type == sf::Event::Closed)
            {
                window.close();
            }
        }

        global_bounds.set_bounds(-view_width / 2 + view_center.x, -view_width / 2 + view_center.y, view_width, view_width);

        window.clear(BACKGROUND_COLOR);
        renderParticles(particles);

        if (track_particle)
            renderParticleInfo(track_particle);

        renderTime(time);
        if (shouldRenderTree)
            renderTree(tree);
        window.display();
    }
}

void Render::findTrackParticle()
{
    sf::Vector2i mouse = sf::Mouse::getPosition(window);
    vector2D mouseLocation = invtransform(vector2D(static_cast<double>(mouse.x), static_cast<double>(mouse.y)));

    query_bounds.set_bounds(mouseLocation.x - 0.05 * view_width, mouseLocation.y - 0.05 * view_width, 0.1 * view_width, 0.1 * view_width);

    query_particles.clear();
    tree->query(query_bounds, query_particles);
    double dist = 1e10;
    for (auto const &particle : query_particles)
    {
        vector2D diff = mouseLocation - particle->position;
        if (diff.norm() < dist)
        {
            dist = diff.norm();
            track_particle = particle;
        }
    }
}