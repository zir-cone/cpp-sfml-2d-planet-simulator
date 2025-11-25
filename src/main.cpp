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
// bodies
// --------------------
struct BodyTypeInfo {
    std::string name;
    std::string description;
};

static BodyTypeInfo genericPlanetType {
    "cold arid nontectonic subearth",
    "non-tectonic rocky planet with low surface temperatures and thin/absent atmosphere"
};

struct Body {
    double mass;
    Vec2 position;
    Vec2 velocity;
    Vec2 acceleration;
    sf::Color color;
    // phys properties
    double density;
    double radius;
    bool fixed = false;
    // metadate
    std::string name;
    BodyTypeInfo typeInfo = genericPlanetType;
    // trail
    std::vector<Vec2> trail;
};

// --------------------
// Global sim parameters
// --------------------

// grav const
// not realistic, for now

const double G = 200.0;
const double SOFTENING = 10.0;

bool trailsEnabled = true;
int maxTrailPoints = 200;
// volume lock
bool lockVolume = false;
// helpers for radius/density
double computeVolumeFromRadius(double r) {
    const double pi = 3.141592653589793;
    return 4.0/3.0*pi*r*r*r;
}

double computeRadiusFromMassDensity(double mass, double density) {
    if (density <= 0.0) density = 1.0;
    double volume = mass / density;
    const double pi = 3.141592653589793;
    double r = std::cbrt((3.0 * volume) / (4.0 * pi));
    return r;
}

void updateBodyAfterMassChange(Body& b) {
    if (!lockVolume) {
        // keep dens and recompute radius
        b.radius = computeRadiusFromMassDensity(b.mass, b.density);
    } else {
        // keep radius, recompute dens
        double volume = computeVolumeFromRadius(b.radius);
        if (volume <= 0.0) volume = 1.0;
        b.density = b.mass / volume;
    }
}

// --------------------
// Compute grav accelerations
// --------------------
void computeGravity(std::vector<Body>& bodies) {
    std::size_t n = bodies.size();
    for (auto& b : bodies) {
        b.acceleration = {0.0, 0.0};
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
// integrator
// this shit is so fucking painful i just want to sleep
// --------------------

void step(std::vector<Body>& bodies, double dt) {
    // update velocities
    for (auto& b : bodies) {
        if (!b.fixed) {
            b.velocity += b.acceleration * dt;
        }
    }
    // update pos
    for (auto& b : bodies) {
        if (!b.fixed) {
            b.position += b.velocity * dt;
        }
    }
    computeGravity(bodies);
    // update trails
    if (trailsEnabled) {
        for (auto& b : bodies) {
            b.trail.push_back(b.position);
            if ((int)b.trail.size() > maxTrailPoints) {
                b.trail.erase(b.trail.begin());
            }
        }
    } else {
        for (auto& b : bodies) b.trail.clear();
    }
}

// --------------------
// interaction modes
// 
// late registration is a classic
// --------------------
enum class InteractionMode {
    AddStill,
    AddMoving,
    AddOrbiting
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

// camera rel body
Vec2 viewCenterWorld{0.0, 0.0};
int viewBodyIndex = -1;

sf::Vector2f worldToScreen(const Vec2& p, double width, double height) {
    double cx = width * 0.5;
    double cy = height * 0.5;
    double x = p.x - viewCenterWorld.x + cx;
    double y = p.y - viewCenterWorld.y + cy;
    return sf::Vector2f(static_cast<float>(x), static_cast<float>(y));
}

// gui panel for sel body
struct BodyPanel {
    bool visible = false;
    bool minimized = false;
    sf::Vector2f pos{20.0f, 80.0f};
    sf::Vector2f size{260.0f, 180.0f};
    bool dragging = false;
    sf::Vector2f dragOffset{0.0f, 0.0f};
};
bool pointInRect(const sf::Vector2f& p, const sf::FloatRect& r) {
    return r.contains(p);
}

int main() {
    const unsigned int WIDTH = 1280;
    const unsigned int HEIGHT = 720;

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "PlanetSim.exe");
    window.setFramerateLimit(120);

    std::vector<Body> bodies;

    // create star
    Body star;
    star.mass = 4000.0;
    star.position = {WIDTH / 2.0, HEIGHT / 2.0 };
    star.velocity = { 0.0, 0.0 };
    star.color = sf::Color::Yellow;
    star.density = 1.0;
    star.radius = 20.0;
    star.fixed = true;
    star.name = "Sol";
    bodies.push_back(star);

    computeGravity(bodies);

    InteractionMode mode = InteractionMode::AddStill;
    int selectedIndex = -1;

    bool isDraggingAddMoving = false;
    Vec2 dragStart{0.0, 0.0};
    sf::Clock clock;
    sf::Font font;
    bool fontLoaded = font.loadFromFile("resources/arial.ttf");
    
    BodyPanel bodyPanel;

    bool renaming = false;
    std::string renameBuffer;

    viewCenterWorld = { WIDTH / 2.0, HEIGHT / 2.0 };

    while (window.isOpen()) {
        sf::Event event;
        while(window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Num1) {
                    mode = InteractionMode::AddStill;
                    std::cout << "Mode: Add Still\n";
                } else if (event.key.code == sf::Keyboard::Num2) {
                    mode = InteractionMode::AddMoving;
                    std::cout << "Mode: Add Moving\n";
                } else if (event.key.code == sf::Keyboard::Num3) {
                    mode = InteractionMode::AddOrbiting;
                    std::cout << "Mode: Add Orbiting\n";
                } else if (event.key.code == sf::Keyboard::Delete) {
                    if (selectedIndex >= 0 && selectedIndex < (int)bodies.size()) {
                        bodies.erase(bodies.begin() + selectedIndex);
                        selectedIndex = -1;
                        bodyPanel.visible = false;
                        computeGravity(bodies);
                    }
                } else if (event.key.code == sf::Keyboard::Up) {
                    if (selectedIndex >= 0) {
                        bodies[selectedIndex].mass *= 1.2;
                        updateBodyAfterMassChange(bodies[selectedIndex]);
                    }
                } else if (event.key.code == sf::Keyboard::Down) {
                    if (selectedIndex >= 0) {
                        bodies[selectedIndex].mass *= 0.8;
                        updateBodyAfterMassChange(bodies[selectedIndex]);
                    }
                } else if (event.key.code == sf::Keyboard::Space) {
                    bodies.clear();
                    bodies.push_back(star);
                    computeGravity(bodies);
                    selectedIndex = -1;
                    bodyPanel.visible = false;
                    viewBodyIndex = -1;
                    viewCenterWorld = { WIDTH / 2.0, HEIGHT / 2.0};
                    for (auto& b : bodies) b.trail.clear();
                } else if (event.key.code == sf::Keyboard::T) {
                    trailsEnabled = !trailsEnabled;
                } else if (event.key.code == sf::Keyboard::LBracket) {
                    if (maxTrailPoints > 10) maxTrailPoints -= 10;
                } else if (event.key.code == sf::Keyboard::RBracket) {
                    maxTrailPoints += 10;
                } else if (event.key.code == sf::Keyboard::V) {
                    lockVolume = !lockVolume;
                } else if (event.key.code == sf::Keyboard::F) {
                    if (selectedIndex >= 0) {
                        if (viewBodyIndex == selectedIndex) {
                            viewBodyIndex = -1;
                        } else {
                            viewBodyIndex = selectedIndex;
                        }
                    }
                } else if (event.key.code == sf::Keyboard::X) {
                    if (selectedIndex >= 0) {
                        bodies[selectedIndex].fixed = !bodies[selectedIndex].fixed;
                    }
                } else if (event.key.code == sf::Keyboard::F2) {
                    if (selectedIndex >= 0) {
                        renaming = true;
                        renameBuffer = bodies[selectedIndex].name;
                    }
                } else if (event.key.code == sf::Keyboard::Escape) {
                    renaming = false;
                }
            }
            // text input for renaming
            else if (event.type == sf::Event::TextEntered && renaming) {
                if (event.text.unicode == '\r' || event.text.unicode == '\n') {
                    // enter: commit
                    if (selectedIndex >= 0) {
                        bodies[selectedIndex].name = renameBuffer;
                    }
                    renaming = false;
                } else if (event.text.unicode == 8) {
                    // backspace
                    if (!renameBuffer.empty()) renameBuffer.pop_back();
                } else if (event.text.unicode >= 32 && event.text.unicode < 127) {
                    renameBuffer.push_back(static_cast<char>(event.text.unicode));
                }
            }
            // mouse press
            else if (event.type == sf::Event::MouseButtonPressed) {
                sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));

                if (bodyPanel.visible) {
                    sf::FloatRect panelRect(bodyPanel.pos, bodyPanel.size);
                    sf::FloatRect titleRect(bodyPanel.pos.x, bodyPanel.pos.y, bodyPanel.size.x, 24.0f);

                    if (event.mouseButton.button == sf::Mouse::Left) {
                        if (titleRect.contains(mouseScreenF)) {
                            float btnSize = 16.0f;
                            sf::FloatRect closeRect(
                                bodyPanel.pos.x + bodyPanel.size.x - btnSize - 4.0f,
                                bodyPanel.pos.y + 4.0f,
                                btnSize, btnSize
                            );
                            sf::FloatRect minRect(
                                bodyPanel.pos.x + bodyPanel.size.x - 2.0f * (btnSize + 4.0f),
                                bodyPanel.pos.y + 4.0f,
                                btnSize, btnSize
                            );

                            if (closeRect.contains(mouseScreenF)) {
                                // deselect on close
                                bodyPanel.visible = false;
                                selectedIndex = -1;
                                continue;
                            } else if (minRect.contains(mouseScreenF)) {
                                bodyPanel.minimized = !bodyPanel.minimized;
                                continue;
                            } else {
                                // dragging
                                bodyPanel.dragging = true;
                                bodyPanel.dragOffset = mouseScreenF - bodyPanel.pos;
                                continue;
                            }
                        } else if (panelRect.contains(mouseScreenF)) {
                            // focus on panel
                            continue;
                        }
                    }
                }

                // panel doesn't want the event or we clicked outside it
                if (event.mouseButton.button == sf::Mouse::Left) {
                    if (mode == InteractionMode::AddStill || mode == InteractionMode::AddMoving || mode == InteractionMode::AddOrbiting) {
                        // Convert mouse to world coordinates based on viewCenterWorld
                        double cx = WIDTH * 0.5;
                        double cy = HEIGHT * 0.5;
                        Vec2 worldPos{
                            mouseScreenF.x - cx + viewCenterWorld.x,
                            mouseScreenF.y - cy + viewCenterWorld.y
                        };

                        if (mode == InteractionMode::AddStill) {
                            Body b;
                            b.mass = 50.0;
                            b.density = 1.0;
                            b.radius = computeRadiusFromMassDensity(b.mass, b.density);
                            b.position = worldPos;
                            b.velocity = {0.0, 0.0};
                            b.color = sf::Color::Cyan;
                            b.name = "Body";
                            bodies.push_back(b);
                            computeGravity(bodies);
                        } else if (mode == InteractionMode::AddMoving) {
                            isDraggingAddMoving = true;
                            dragStart = worldPos;
                        } else if (mode == InteractionMode::AddOrbiting) {
                            if (selectedIndex >= 0) {
                                Body b;
                                b.mass = 50.0;
                                b.density = 1.0;
                                b.radius = computeRadiusFromMassDensity(b.mass, b.density);
                                b.position = worldPos;
                                b.color = sf::Color::Green;
                                b.name = "Orbiting body";

                                Body& parent = bodies[selectedIndex];
                                Vec2 r = b.position - parent.position;
                                double dist = length(r);
                                if (dist < 1.0) dist = 1.0;
                                double vmag = std::sqrt(G * parent.mass / dist);

                                // Perpendicular direction
                                Vec2 perp{-r.y / dist, r.x / dist};
                                Vec2 vOrb = perp * vmag;

                                b.velocity = parent.velocity + vOrb;

                                bodies.push_back(b);
                                computeGravity(bodies);
                            }
                        }
                    }
                } else if (event.mouseButton.button == sf::Mouse::Right) {
                    // select nearest body
                    double cx = WIDTH * 0.5;
                    double cy = HEIGHT * 0.5;
                    Vec2 worldPos{
                        mouseScreenF.x - cx + viewCenterWorld.x,
                        mouseScreenF.y - cy + viewCenterWorld.y
                    };
                    sf::Vector2f worldF(
                        static_cast<float>(worldPos.x),
                        static_cast<float>(worldPos.y)
                    );
                    int idx = findNearestBody(bodies, worldF, 25.0f);
                    selectedIndex = idx;
                    if (idx >= 0) {
                        bodyPanel.visible = true;
                        bodyPanel.minimized = false;
                    }
                }
            }

            // mouse release
            else if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    if (bodyPanel.dragging) {
                        bodyPanel.dragging = false;
                    }
                    if (mode == InteractionMode::AddMoving && isDraggingAddMoving) {
                        isDraggingAddMoving = false;

                        sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                        double cx = WIDTH * 0.5;
                        double cy = HEIGHT * 0.5;
                        Vec2 worldEnd{
                            mouseScreenF.x - cx + viewCenterWorld.x,
                            mouseScreenF.y - cy + viewCenterWorld.y
                        };

                        Vec2 vel = worldEnd - dragStart;
                        vel *= 0.5;

                        Body b;
                        b.mass = 50.0;
                        b.density = 1.0;
                        b.radius = computeRadiusFromMassDensity(b.mass, b.density);
                        b.position = dragStart;
                        b.velocity = vel;
                        b.color = sf::Color::Green;
                        b.name = "Body";
                        bodies.push_back(b);
                        computeGravity(bodies);
                    }
                }
            }

            // --- Mouse move (for dragging panel) ---
            else if (event.type == sf::Event::MouseMoved) {
                if (bodyPanel.dragging) {
                    sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                    bodyPanel.pos = mouseScreenF - bodyPanel.dragOffset;
                }
            }
        } // end pollEvent loop

        // --- Update physics ---
        float dtSeconds = clock.restart().asSeconds();
        if (dtSeconds > 0.05f) dtSeconds = 0.05f;
        double dt = static_cast<double>(dtSeconds);

        if (!bodies.empty()) {
            step(bodies, dt);
        }

        // Update view center
        if (viewBodyIndex >= 0 && viewBodyIndex < (int)bodies.size()) {
            viewCenterWorld = bodies[viewBodyIndex].position;
        }

        // --- Render ---
        window.clear(sf::Color(10, 10, 20));

        // Draw trails
        if (trailsEnabled) {
            for (const auto& b : bodies) {
                if (b.trail.size() < 2) continue;
                std::vector<sf::Vertex> lineVertices;
                lineVertices.reserve(b.trail.size());
                for (const auto& p : b.trail) {
                    sf::Vector2f s = worldToScreen(p, WIDTH, HEIGHT);
                    lineVertices.emplace_back(s, sf::Color(150, 150, 150));
                }
                window.draw(lineVertices.data(), lineVertices.size(), sf::LineStrip);
            }
        }

        // Draw bodies
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            const Body& b = bodies[i];
            sf::Vector2f screenPos = worldToScreen(b.position, WIDTH, HEIGHT);
            float radius = static_cast<float>(b.radius);

            sf::CircleShape circle(radius);
            circle.setOrigin(radius, radius);
            circle.setPosition(screenPos);
            circle.setFillColor(b.color);

            if ((int)i == selectedIndex) {
                circle.setOutlineThickness(2.0f);
                circle.setOutlineColor(sf::Color::Red);
            } else if (b.fixed) {
                circle.setOutlineThickness(1.5f);
                circle.setOutlineColor(sf::Color(200, 200, 50));
            }

            window.draw(circle);
        }

        // HUD
        if (fontLoaded) {
            sf::Text text;
            text.setFont(font);
            text.setCharacterSize(14);
            text.setFillColor(sf::Color::White);
            text.setPosition(10.0f, 10.0f);

            std::string modeStr;
            if (mode == InteractionMode::AddStill) modeStr = "Add Still (1)";
            else if (mode == InteractionMode::AddMoving) modeStr = "Add Moving (2)";
            else modeStr = "Add Orbiting (3)";

            std::string selStr = "None";
            if (selectedIndex >= 0 && selectedIndex < (int)bodies.size()) {
                selStr = bodies[selectedIndex].name;
            }

            std::string hud =
                "M1: add body   M2: select body\n"
                "1: add still   2: add moving   3: add orbiting\n"
                "Del: delete selected   Up/Down: change mass\n"
                "T: toggle trails  [ ]: trail length\n"
                "V: toggle volume lock   X: toggle fixed\n"
                "F: follow selected   F2: rename   Space: reset\n" +
                std::string("Mode: ") + modeStr +
                "   Selected: " + selStr +
                (lockVolume ? "   [Volume LOCKED]" : "   [Volume UNLOCKED]");

            text.setString(hud);
            window.draw(text);
        }

        // Body panel
        if (bodyPanel.visible && selectedIndex >= 0 && selectedIndex < (int)bodies.size() && fontLoaded) {
            Body& b = bodies[selectedIndex];

            sf::RectangleShape panel;
            panel.setPosition(bodyPanel.pos);
            panel.setSize(bodyPanel.size);
            panel.setFillColor(sf::Color(20, 20, 35, 230));
            panel.setOutlineThickness(1.0f);
            panel.setOutlineColor(sf::Color::White);
            window.draw(panel);

            // Title bar
            float titleH = 24.0f;
            sf::RectangleShape titleBar(sf::Vector2f(bodyPanel.size.x, titleH));
            titleBar.setPosition(bodyPanel.pos);
            titleBar.setFillColor(sf::Color(40, 40, 60));
            window.draw(titleBar);

            sf::Text title;
            title.setFont(font);
            title.setCharacterSize(14);
            title.setFillColor(sf::Color::White);
            title.setPosition(bodyPanel.pos.x + 6.0f, bodyPanel.pos.y + 4.0f);
            title.setString("Body Editor");
            window.draw(title);

            // Minimize & close buttons
            float btnSize = 16.0f;
            sf::RectangleShape closeBtn(sf::Vector2f(btnSize, btnSize));
            closeBtn.setPosition(bodyPanel.pos.x + bodyPanel.size.x - btnSize - 4.0f, bodyPanel.pos.y + 4.0f);
            closeBtn.setFillColor(sf::Color(120, 40, 40));
            window.draw(closeBtn);

            sf::RectangleShape minBtn(sf::Vector2f(btnSize, btnSize));
            minBtn.setPosition(bodyPanel.pos.x + bodyPanel.size.x - 2.0f * (btnSize + 4.0f), bodyPanel.pos.y + 4.0f);
            minBtn.setFillColor(sf::Color(80, 80, 80));
            window.draw(minBtn);

            sf::Text xText;
            xText.setFont(font);
            xText.setCharacterSize(12);
            xText.setFillColor(sf::Color::White);
            xText.setPosition(closeBtn.getPosition().x + 4.0f, closeBtn.getPosition().y - 2.0f);
            xText.setString("x");
            window.draw(xText);

            sf::Text mText;
            mText.setFont(font);
            mText.setCharacterSize(12);
            mText.setFillColor(sf::Color::White);
            mText.setPosition(minBtn.getPosition().x + 4.0f, minBtn.getPosition().y - 3.0f);
            mText.setString("-");
            window.draw(mText);

            if (!bodyPanel.minimized) {
                float y = bodyPanel.pos.y + titleH + 6.0f;
                float x = bodyPanel.pos.x + 8.0f;

                sf::Text line;
                line.setFont(font);
                line.setCharacterSize(13);
                line.setFillColor(sf::Color::White);

                // Name
                line.setPosition(x, y);
                if (renaming) {
                    line.setString("Name: " + renameBuffer + "_");
                } else {
                    line.setString("Name: " + b.name + "  (F2 to rename)");
                }
                window.draw(line);
                y += 20.0f;

                // Type
                line.setPosition(x, y);
                line.setString("Type: " + b.typeInfo.name);
                window.draw(line);
                y += 20.0f;

                // Mass
                line.setPosition(x, y);
                line.setString("Mass: " + std::to_string(b.mass));
                window.draw(line);
                y += 20.0f;

                // Density
                line.setPosition(x, y);
                line.setString("Density: " + std::to_string(b.density));
                window.draw(line);
                y += 20.0f;

                // Radius
                line.setPosition(x, y);
                line.setString("Radius: " + std::to_string(b.radius));
                window.draw(line);
                y += 20.0f;

                // Flags
                line.setPosition(x, y);
                std::string flags = "Fixed: ";
                flags += (b.fixed ? "Yes" : "No");
                flags += "   Volume lock: ";
                flags += (lockVolume ? "On" : "Off");
                line.setString(flags);
                window.draw(line);
                y += 20.0f;
            }
        }

        window.display();
    }

    return 0;
}