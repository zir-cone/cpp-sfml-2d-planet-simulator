#include <SFML/Graphics.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

// --------------------
// Basic 2d vector
// --------------------
struct Vec2 {
    double x, y;
    Vec2 operator+(const Vec2& other) const { return {x + other.x, y + other.y}; }
    Vec2 operator-(const Vec2& other) const { return {x - other.x, y - other.y}; }
    Vec2 operator*(double s) const { return {x*s, y*s}; }

    Vec2 operator+=(const Vec2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }
    Vec2 operator *=(double s) {
        x *= s;
        y *= s;
        return *this;
    }
};

double length (const Vec2& v) {
    return std::sqrt(v.x * v.x + v.y * v.y);
}

// --------------------
// Body in the simulation
// --------------------
struct Body {
    double mass;
    Vec2 position;
    Vec2 velocity;
    Vec2 acceleration;
    sf::Color color;
};

// --------------------
// Global sim parameters
// --------------------

// "Game scale" grav const
// Not realistic, for now

const double G = 200.0;
const double SOFTENING = 10.0;

// --------------------
// Compute grav accelerations
// --------------------
void computeGravity(std::vector<Body>& bodies) {
    std::size_t n = bodies.size();
    for (std::size_t i = 0; i < n; ++i) {
        bodies[i].acceleration = {0.0, 0.0};
    }

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            Vec2 r = bodies[j].position - bodies[i].position;
            double dist = length(r) + 1e-6; //avoid div by 0
            double invDist3 = 1.0 / std::pow(dist * dist + SOFTENING * SOFTENING, 1.5);

            Vec2 dir = r;
            double factor = G *invDist3;

            Vec2 accel_i = dir * (factor * bodies[j].mass);
            Vec2 accel_j = dir * (-factor *bodies[i].mass);

            bodies[i].acceleration += accel_i;
            bodies[j].acceleration += accel_j;
        }
    }    
}

// --------------------
// integrator (simple sympletic)
// this shit is so fucking painful i just want to sleep
// --------------------
void step(std::vector<Body>& bodies, double dt) {
    for (auto& b : bodies) {
        b.velocity += b.acceleration * dt;
    }
    for (auto& b : bodies) {
        b.position += b.velocity * dt;
    }

    computeGravity(bodies);
}

// --------------------
// interaction modes
// 
// late registration is a classic
// --------------------
enum class InteractionMode {
    AddStill,
    AddMoving
};
// find nearest body to a point on the screen
int findNearestBody(const std::vector<Body>& bodies, const sf::Vector2f& point, float maxDistPixels) {
    int bestIndex = -1;
    double bestDist2 = static_cast<double>(maxDistPixels) * maxDistPixels;

    for (std::size_t i = 0; i < bodies.size(); ++i) {
        double dx = bodies[i].position.x - point.x;
        double dy = bodies [i].position.y - point.y;
        double d2 = dx*dx + dy*dy;
        if (d2 < bestDist2) {
            bestDist2 = d2;
            bestIndex = static_cast<int>(i);
        }
    }
    return bestIndex;
}

// convert a mass into a visible radius for drawing
float radiusFromMass(double mass) {
    float r = static_cast<float>(std::cbrt(mass)); // cube root is sexy
    if (r < 3.0f) r = 3.0f;
    if (r > 40.0f) r = 40.0f;
    return r;
}
// FUCK YOU IM GOING TO BED
// i couldnt sleep :wilted_rose:

int main() {
    const unsigned int WIDTH = 1280;
    const unsigned int HEIGHT = 720;

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "PlanetSim.exe");
    window.setFramerateLimit(120);

    std::vector<Body> bodies;

    // add a central star, but if a star was just a dingleberry
    Body star;
    star.mass = 2000.0;
    star.position = {WIDTH/2.0, HEIGHT/2.0};
    star.velocity = {0.0, 0.0};
    star.color = sf::Color::Yellow;
    bodies.push_back(star);
    
    computeGravity(bodies);

    InteractionMode mode = InteractionMode::AddStill;
    int selectedIndex = -1;

    bool isDragging = false;
    Vec2 dragStart {0.0, 0.0};

    sf::Clock clock;

    // font that ill handle ladder
    sf::Font font;
    bool fontLoaded = font.loadFromFile("resources/arial.ttf");

    while (window.isOpen()) {
        // handle events
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Num1) {
                    mode = InteractionMode::AddStill;
                    std::cout << "Mode: Add Still\n";
                }
                else if (event.key.code == sf::Keyboard::Num2) {
                    mode = InteractionMode::AddMoving;
                    std::cout << "Mode: Add Moving\n";
                }
                else if (event.key.code == sf::Keyboard::Delete) {
                    if (selectedIndex >= 0 && selectedIndex < static_cast<int>(bodies.size())) {
                        bodies.erase(bodies.begin() + selectedIndex);
                        selectedIndex = -1;
                        computeGravity(bodies);
                        std::cout << "Delete selected body.\n";
                    }
                }
                else if (event.key.code == sf::Keyboard::Up) { //UP!!
                    if (selectedIndex >= 0) {
                        bodies[selectedIndex].mass *= 1.2;
                        std::cout << "Increased mass to " << bodies[selectedIndex].mass << "\n";

                    }
                }
                else if (event.key.code == sf::Keyboard::Down) {
                    if (selectedIndex >= 0) {
                        bodies[selectedIndex].mass *= 0.8;
                        std::cout << "Decreased mass to " << bodies[selectedIndex].mass << "\n";
                    }
                }
                else if (event.key.code == sf::Keyboard::Space) {
                    bodies.clear();
                    bodies.push_back(star);
                    computeGravity(bodies);
                    selectedIndex = -1;
                    std::cout << "Reset system.\n";
                }
            }
            else if (event.type == sf::Event::MouseButtonPressed) {
                sf::Vector2 mousePosF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                if (event.mouseButton.button == sf::Mouse::Left) {
                    if (mode == InteractionMode::AddStill) {
                        Body b;
                        b.mass = 50.0;
                        b.position = { mousePosF.x, mousePosF.y };
                        b.velocity = { 0.0, 0.0 };
                        b.color = sf::Color::Cyan;
                        bodies.push_back(b);
                        computeGravity(bodies);
                        std::cout << "Add still body at (" << b.position.x << ", " << b.position.y << ")\n";
                    }
                    else if (mode == InteractionMode::AddMoving) {
                        isDragging = true;
                        dragStart = { mousePosF.x, mousePosF.y };
                    }
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    // select nearest body
                    int idx = findNearestBody(bodies, mousePosF, 30.0f);
                    selectedIndex = idx;
                    if (idx >= 0) {
                        std::cout << "Selected body #" << idx
                            << "pos=(" << bodies[idx].position.x << "," << bodies[idx].position.y << ")"
                            << "mass=" << bodies[idx].mass << "\n";    
                    } else {
                        std::cout << "No body near that click!\n";
                    }
                }
            }
            else if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left &&
                mode == InteractionMode::AddMoving && isDragging) {
                    isDragging = false;
                    sf::Vector2f mousePosF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                    Vec2 dragEnd { mousePosF.x, mousePosF.y };

                    Vec2 vel = dragEnd - dragStart;
                    vel *= 0.5;
                    Body b;
                    b.mass = 50.0;
                    b.position = dragStart;
                    b.velocity = vel;
                    b.color = sf::Color::Green;
                    bodies.push_back(b);
                    computeGravity(bodies);
                    std::cout << "Added moving body at (" << b.position.x << "," << b.position.y
                        << ") vel=(" << b.velocity.x << "," << b.velocity.y << ")\n";
                }
            }
        }
        // --------------------
        // update physics
        // --------------------
        float dtSeconds = clock.restart().asSeconds();
        if (dtSeconds > 0.05f) dtSeconds = 0.05f;
        double dt = static_cast<double>(dtSeconds);

        if (!bodies.empty()) {
            step(bodies, dt);
        }
        // --------------------
        // render
        // --------------------
        window.clear(sf::Color(10, 10, 20));
        // draw vel vector
        if (isDragging && mode == InteractionMode::AddMoving) {
            sf::Vector2f mousePosF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
            sf::Vertex line[] = {
                sf::Vertex(sf::Vector2f(dragStart.x, dragStart.y), sf::Color::White),
                sf::Vertex(mousePosF, sf::Color::White)
            };
            window.draw(line, 2, sf::Lines);
        }
        // --------------------
        // draw bodies
        // --------------------
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            const Body& b = bodies[i];
            float radius = radiusFromMass(b.mass);
            sf::CircleShape circle(radius);
            circle.setOrigin(radius, radius);
            circle.setPosition(static_cast<float>(b.position.x),
            static_cast<float>(b.position.y));
            circle.setFillColor(b.color);
            if (static_cast<int>(i) == selectedIndex) {
                circle.setOutlineThickness(2.0f);
                circle.setOutlineColor(sf::Color::Red);
            }
            window.draw(circle);
        }
        // gui txt
        if (fontLoaded) {
            sf::Text text;
            text.setFont(font);
            text.setCharacterSize(14);
            text.setFillColor(sf::Color::White);
            text.setPosition(10.0f, 10.0f);
            
            std::string modeStr =
            (mode == InteractionMode::AddStill) ? "Add Still (1)" : "Add Moving (2)";
            std::string selStr = "None";
            if (selectedIndex >= 0 && selectedIndex < static_cast<int>(bodies.size())) {
                const Body& b = bodies[selectedIndex];
                selStr = "idx=" + std::to_string(selectedIndex) +
                " m=" + std::to_string(b.mass);
            }
            text.setString(
                "M1: add body   M2: select body\n"
                "1: add still   2 add moving    Del: delete selected\n"
                "Up/Down: change mass of selected   Space: reset\n"
                "Mode: " + modeStr + "  Selected: " +selStr
            );
            window.draw(text);
        }
        window.display();
    }
    return 0;
}

// all ts 11/23/25 9:34 P.M --> 11/24/25 2:34 A.M.
// im SO tireddd