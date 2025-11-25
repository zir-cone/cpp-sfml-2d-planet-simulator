#include <SFML/Graphics.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <cstdio>   // snprintf

// --------------------
// Basic 2d vector
// --------------------
struct Vec2 {
    double x, y;

    Vec2 operator+(const Vec2& other) const { return {x + other.x, y + other.y}; }
    Vec2 operator-(const Vec2& other) const { return {x - other.x, y - other.y}; }
    Vec2 operator*(double s)        const { return {x * s, y * s}; }

    Vec2& operator+=(const Vec2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }
    Vec2& operator*=(double s) {
        x *= s;
        y *= s;
        return *this;
    }
};

double length(const Vec2& v) {
    return std::sqrt(v.x * v.x + v.y * v.y);
}

// --------------------
// Body types
// --------------------
struct BodyTypeInfo {
    std::string name;
    std::string description;
};

static BodyTypeInfo genericPlanetType{
    "cold arid nontectonic subearth",
    "non-tectonic rocky planet with low surface temperatures and thin/absent atmosphere"
};

// --------------------
// Bodies
// --------------------
struct Body {
    double mass;         // kg
    Vec2   position;     // meters
    Vec2   velocity;     // m/s
    Vec2   acceleration; // m/s^2
    sf::Color color;

    double density;
    double radius;       // pixels
    bool   fixed = false;
    bool   ghost = false; // ghost: visible but no gravity / motion

    std::string  name;
    BodyTypeInfo typeInfo = genericPlanetType;

    std::vector<Vec2> trail; // positions in meters
};

// --------------------
// Global sim parameters
// --------------------

// Physical constants
const double G_PHYS = 6.67430e-11;          // m^3 kg^-1 s^-2
const double YEAR_SECONDS = 365.25 * 24.0 * 3600.0;

const double M_SUN   = 1.98847e30;         // kg
const double M_EARTH = 5.97219e24;         // kg

// Simulation scale (mutable)
double metersPerPixel = 1.0e9;             // 1 px = 1e9 m
double timeScale      = 86400.0;           // 1 real s = 1 simulated day

// Softening in meters
const double SOFTENING_METERS = 1.0e9;

bool trailsEnabled = true;
int  maxTrailPoints = 200;
bool lockVolume     = false;               // volume lock for mass editing

// helpers for radius/density (just internal consistency)
double computeVolumeFromRadius(double r) {
    const double pi = 3.141592653589793;
    return 4.0 / 3.0 * pi * r * r * r;
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
        b.radius = computeRadiusFromMassDensity(b.mass, b.density);
    } else {
        double volume = computeVolumeFromRadius(b.radius);
        if (volume <= 0.0) volume = 1.0;
        b.density = b.mass / volume;
    }
}

// simple scientific notation formatter
std::string formatSci(double v, int prec = 3) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*e", prec, v);
    return std::string(buf);
}

// nice distance formatter
std::string formatDistance(double meters) {
    char buf[64];
    if (meters >= 1.0e12)      std::snprintf(buf, sizeof(buf), "%.3f Tm", meters / 1.0e12);
    else if (meters >= 1.0e9)  std::snprintf(buf, sizeof(buf), "%.3f Gm", meters / 1.0e9);
    else if (meters >= 1.0e6)  std::snprintf(buf, sizeof(buf), "%.3f Mm", meters / 1.0e6);
    else if (meters >= 1.0e3)  std::snprintf(buf, sizeof(buf), "%.3f km", meters / 1.0e3);
    else                       std::snprintf(buf, sizeof(buf), "%.0f m", meters);
    return std::string(buf);
}

// --------------------
// Gravity (Newtonian)
// --------------------
void computeGravity(std::vector<Body>& bodies) {
    std::size_t n = bodies.size();

    for (auto& b : bodies) {
        b.acceleration = {0.0, 0.0};
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (bodies[i].ghost) continue; // ghosts don't feel gravity
        for (std::size_t j = i + 1; j < n; ++j) {
            if (bodies[j].ghost) continue; // ghosts don't exert gravity

            Vec2 r = bodies[j].position - bodies[i].position; // meters
            double dist2 = r.x * r.x + r.y * r.y;
            double softened2 = dist2 + SOFTENING_METERS * SOFTENING_METERS;
            double softened = std::sqrt(softened2);
            if (softened < 1.0) softened = 1.0;

            double invDist3 = 1.0 / (softened2 * softened);

            Vec2 accel_i = r * (G_PHYS * bodies[j].mass * invDist3);
            Vec2 accel_j = r * (-G_PHYS * bodies[i].mass * invDist3);

            bodies[i].acceleration += accel_i;
            bodies[j].acceleration += accel_j;
        }
    }
}

// --------------------
// integrator
// --------------------
void step(std::vector<Body>& bodies, double dtRealSeconds) {
    double dt = dtRealSeconds * timeScale;  // seconds of simulated time

    for (auto& b : bodies) {
        if (!b.fixed && !b.ghost) {
            b.velocity += b.acceleration * dt;
        }
    }

    for (auto& b : bodies) {
        if (!b.fixed && !b.ghost) {
            b.position += b.velocity * dt;
        }
    }

    computeGravity(bodies);

    if (trailsEnabled) {
        for (auto& b : bodies) {
            if (!b.ghost) {
                b.trail.push_back(b.position);
                if ((int)b.trail.size() > maxTrailPoints) {
                    b.trail.erase(b.trail.begin());
                }
            }
        }
    } else {
        for (auto& b : bodies) b.trail.clear();
    }
}

// --------------------
// interaction modes
// --------------------
enum class InteractionMode {
    AddStill,
    AddMoving,
    AddOrbitingPlace, // click to place ghost
    OrbitEdit         // ghost orbit edit + overlay
};

// find nearest body to a world-space point (meters)
int findNearestBody(const std::vector<Body>& bodies, const Vec2& pointMeters, double maxDistMeters) {
    int bestIndex = -1;
    double bestDist2 = maxDistMeters * maxDistMeters;

    for (std::size_t i = 0; i < bodies.size(); ++i) {
        double dx = bodies[i].position.x - pointMeters.x;
        double dy = bodies[i].position.y - pointMeters.y;
        double d2 = dx * dx + dy * dy;
        if (d2 < bestDist2) {
            bestDist2 = d2;
            bestIndex = static_cast<int>(i);
        }
    }
    return bestIndex;
}

// choose parent body by maximum GM/r^2 at a point
int findStrongestInfluenceParent(const std::vector<Body>& bodies, const Vec2& posMeters) {
    int bestIndex = -1;
    double bestAccel = 0.0;

    for (std::size_t i = 0; i < bodies.size(); ++i) {
        if (bodies[i].ghost) continue;
        Vec2 r = bodies[i].position - posMeters;
        double dist2 = r.x * r.x + r.y * r.y;
        if (dist2 <= 0.0) continue;

        double accel = G_PHYS * bodies[i].mass / dist2;
        if (accel > bestAccel) {
            bestAccel = accel;
            bestIndex = static_cast<int>(i);
        }
    }
    return bestIndex;
}

// camera rel body
Vec2 viewCenterWorld{0.0, 0.0};
int  viewBodyIndex = -1;

// world (meters) -> screen (pixels)
sf::Vector2f worldToScreen(const Vec2& p, double width, double height) {
    double cx = width * 0.5;
    double cy = height * 0.5;

    double dx = (p.x - viewCenterWorld.x) / metersPerPixel;
    double dy = (p.y - viewCenterWorld.y) / metersPerPixel;

    double x = cx + dx;
    double y = cy + dy;

    return sf::Vector2f(static_cast<float>(x), static_cast<float>(y));
}

// gui panel for sel body
struct BodyPanel {
    bool visible   = false;
    bool minimized = false;
    sf::Vector2f pos{20.0f, 80.0f};
    sf::Vector2f size{260.0f, 180.0f};
    bool dragging  = false;
    sf::Vector2f dragOffset{0.0f, 0.0f};
};

// --------------------
// Orbit overlay for OrbitEdit
// --------------------
void drawOrbitOverlay(
    sf::RenderWindow& window,
    const Body& body,
    const Body& parent,
    double width,
    double height,
    const sf::Font& font
) {
    // Barycenter (2-body)
    double Mtot = body.mass + parent.mass;
    Vec2 bary{
        (body.position.x * body.mass + parent.position.x * parent.mass) / Mtot,
        (body.position.y * body.mass + parent.position.y * parent.mass) / Mtot
    };

    Vec2 relBody = body.position - bary;
    Vec2 relParent = parent.position - bary;
    double a = length(relBody);
    if (a <= 0.0) return;

    // approximated as circle radius a
    const int SEGMENTS = 256;
    std::vector<sf::Vertex> orbitLine;
    orbitLine.reserve(SEGMENTS + 1);
    for (int i = 0; i <= SEGMENTS; ++i) {
        double t = (double)i / SEGMENTS;
        double ang = t * 2.0 * 3.141592653589793;
        Vec2 pWorld{ bary.x + a * std::cos(ang), bary.y + a * std::sin(ang) };
        sf::Vector2f s = worldToScreen(pWorld, width, height);
        orbitLine.emplace_back(s, sf::Color(0, 150, 255)); // blue
    }
    window.draw(orbitLine.data(), orbitLine.size(), sf::LineStrip);

    // barycenter point
    sf::Vector2f baryScreen = worldToScreen(bary, width, height);
    sf::CircleShape baryPt(4.0f);
    baryPt.setOrigin(4.0f, 4.0f);
    baryPt.setPosition(baryScreen);
    baryPt.setFillColor(sf::Color::Red);
    window.draw(baryPt);

    // Pe direction: from barycenter towards parent
    double distParent = length(relParent);
    if (distParent <= 0.0) return;
    Vec2 peDir{ relParent.x / distParent, relParent.y / distParent };

    Vec2 peWorld{ bary.x + a * peDir.x, bary.y + a * peDir.y };
    Vec2 apWorld{ bary.x - a * peDir.x, bary.y - a * peDir.y };

    sf::Vector2f peScreen = worldToScreen(peWorld, width, height);
    sf::Vector2f apScreen = worldToScreen(apWorld, width, height);

    // dashed lines
    auto drawDashed = [&](const sf::Vector2f& A, const sf::Vector2f& B) {
        const float dashLen = 8.0f;
        const float gapLen  = 4.0f;
        sf::Vector2f d = B - A;
        float len = std::sqrt(d.x*d.x + d.y*d.y);
        if (len <= 0.0f) return;
        sf::Vector2f dir = d / len;
        float pos = 0.0f;
        while (pos < len) {
            float seg = std::min(dashLen, len - pos);
            sf::Vector2f s = A + dir * pos;
            sf::Vector2f e = A + dir * (pos + seg);
            sf::Vertex line[] = {
                sf::Vertex(s, sf::Color::White),
                sf::Vertex(e, sf::Color::White)
            };
            window.draw(line, 2, sf::Lines);
            pos += dashLen + gapLen;
        }
    };

    drawDashed(baryScreen, peScreen);
    drawDashed(baryScreen, apScreen);

    // Pe / Ap markers
    sf::CircleShape pePt(4.0f);
    pePt.setOrigin(4.0f, 4.0f);
    pePt.setPosition(peScreen);
    pePt.setFillColor(sf::Color(255, 100, 200));
    window.draw(pePt);

    sf::CircleShape apPt(4.0f);
    apPt.setOrigin(4.0f, 4.0f);
    apPt.setPosition(apScreen);
    apPt.setFillColor(sf::Color(255, 100, 200));
    window.draw(apPt);

    if (font.getInfo().family.empty()) return;

    sf::Text text;
    text.setFont(font);
    text.setCharacterSize(12);

    text.setFillColor(sf::Color(255, 100, 200));
    text.setString("Pe");
    text.setPosition(peScreen.x + 6.0f, peScreen.y - 4.0f);
    window.draw(text);

    text.setString("Ap");
    text.setPosition(apScreen.x + 6.0f, apScreen.y - 4.0f);
    window.draw(text);

    // distances (circular, so r(Pe)=r(Ap)=a)
    text.setFillColor(sf::Color::White);
    text.setString("r(Pe) = " + formatDistance(a));
    text.setPosition((baryScreen.x + peScreen.x) * 0.5f + 4.0f,
                     (baryScreen.y + peScreen.y) * 0.5f);
    window.draw(text);

    text.setString("r(Ap) = " + formatDistance(a));
    text.setPosition((baryScreen.x + apScreen.x) * 0.5f + 4.0f,
                     (baryScreen.y + apScreen.y) * 0.5f + 14.0f);
    window.draw(text);

    // true anomaly ν (angle from Pe)
    Vec2 relFromBary = body.position - bary;
    double angBody = std::atan2(relFromBary.y, relFromBary.x);
    double angPe   = std::atan2(peDir.y, peDir.x);
    double nu = angBody - angPe;
    while (nu < 0.0) nu += 2.0 * 3.141592653589793;
    while (nu >= 2.0 * 3.141592653589793) nu -= 2.0 * 3.141592653589793;
    double nuDeg = nu * 180.0 / 3.141592653589793;

    char buf[64];
    std::snprintf(buf, sizeof(buf), "ν = %.2f°", nuDeg);
    text.setFillColor(sf::Color::Green);
    text.setString(buf);
    sf::Vector2f bodyScreen = worldToScreen(body.position, width, height);
    text.setPosition(bodyScreen.x + 8.0f, bodyScreen.y - 4.0f);
    window.draw(text);

    // orbital period (Kepler 3rd, 2-body)
    double T = 2.0 * 3.141592653589793 * std::sqrt(a * a * a / (G_PHYS * (body.mass + parent.mass)));
    double T_years = T / YEAR_SECONDS;
    std::snprintf(buf, sizeof(buf), "P = %.3f yr", T_years);
    text.setFillColor(sf::Color(0, 180, 255));
    text.setString(buf);
    sf::Vector2f labelPos = worldToScreen({bary.x + a, bary.y}, width, height);
    text.setPosition(labelPos.x + 6.0f, labelPos.y + 4.0f);
    window.draw(text);
}

// scale / time formatting
std::string formatScale() {
    return formatDistance(metersPerPixel) + " / px";
}
std::string formatTimeScaleStr() {
    double daysPerSec = timeScale / 86400.0;
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.3f days/s", daysPerSec);
    return std::string(buf);
}

int main() {
    const unsigned int WIDTH  = 1280;
    const unsigned int HEIGHT = 720;

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "PlanetSim.exe");
    window.setFramerateLimit(120);

    std::vector<Body> bodies;

    // create star at origin in meters
    Body star;
    star.mass     = M_SUN;
    star.position = {0.0, 0.0};
    star.velocity = {0.0, 0.0};
    star.color    = sf::Color::Yellow;
    star.density  = 1.0;
    star.radius   = 20.0; // pixels
    star.fixed    = true;
    star.name     = "Sol";
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

    // orbit edit state
    int orbitEditBodyIndex   = -1;
    int orbitEditParentIndex = -1;

    // view centered on star
    viewCenterWorld = star.position;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            // -------------------- Events --------------------
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed) {
                // global hotkeys
                if (event.key.code == sf::Keyboard::Num1) {
                    mode = InteractionMode::AddStill;
                } else if (event.key.code == sf::Keyboard::Num2) {
                    mode = InteractionMode::AddMoving;
                } else if (event.key.code == sf::Keyboard::Num3) {
                    mode = InteractionMode::AddOrbitingPlace;
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
                } else if (event.key.code == sf::Keyboard::Space && !renaming) {
                    // reset system
                    bodies.clear();
                    star.position = {0.0, 0.0};
                    star.velocity = {0.0, 0.0};
                    star.fixed    = true;
                    bodies.push_back(star);
                    computeGravity(bodies);
                    selectedIndex = -1;
                    bodyPanel.visible = false;
                    viewBodyIndex = -1;
                    viewCenterWorld = star.position;
                    for (auto& b : bodies) b.trail.clear();
                    orbitEditBodyIndex = orbitEditParentIndex = -1;
                    mode = InteractionMode::AddStill;
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
                    // cancel rename OR orbit edit
                    if (renaming) renaming = false;
                    if (mode == InteractionMode::OrbitEdit &&
                        orbitEditBodyIndex >= 0 &&
                        orbitEditBodyIndex < (int)bodies.size()) {
                        bodies.erase(bodies.begin() + orbitEditBodyIndex);
                        selectedIndex = -1;
                        bodyPanel.visible = false;
                        orbitEditBodyIndex = orbitEditParentIndex = -1;
                        mode = InteractionMode::AddOrbitingPlace;
                        computeGravity(bodies);
                    }
                }
                // zoom + time scale adjustments
                else if (event.key.code == sf::Keyboard::Z) { // zoom in
                    metersPerPixel *= 0.5;
                    if (metersPerPixel < 1.0e7) metersPerPixel = 1.0e7;
                } else if (event.key.code == sf::Keyboard::X) { // zoom out
                    metersPerPixel *= 2.0;
                    if (metersPerPixel > 1.0e12) metersPerPixel = 1.0e12;
                } else if (event.key.code == sf::Keyboard::Comma) { // slower time
                    timeScale *= 0.5;
                    if (timeScale < 3600.0) timeScale = 3600.0; // >= 1 hr/s
                } else if (event.key.code == sf::Keyboard::Period) { // faster time
                    timeScale *= 2.0;
                    if (timeScale > 1.0e7) timeScale = 1.0e7;
                }

                // OrbitEdit commit
                if (mode == InteractionMode::OrbitEdit &&
                    (event.key.code == sf::Keyboard::Enter || event.key.code == sf::Keyboard::Return)) {
                    if (orbitEditBodyIndex >= 0 && orbitEditBodyIndex < (int)bodies.size()) {
                        bodies[orbitEditBodyIndex].ghost = false;
                        computeGravity(bodies);
                    }
                    mode = InteractionMode::AddOrbitingPlace;
                    orbitEditBodyIndex = orbitEditParentIndex = -1;
                }
            }
            // text input for renaming
            else if (event.type == sf::Event::TextEntered && renaming) {
                if (event.text.unicode == '\r' || event.text.unicode == '\n') {
                    if (selectedIndex >= 0) {
                        bodies[selectedIndex].name = renameBuffer;
                    }
                    renaming = false;
                } else if (event.text.unicode == 8) {
                    if (!renameBuffer.empty()) renameBuffer.pop_back();
                } else if (event.text.unicode >= 32 && event.text.unicode < 127) {
                    renameBuffer.push_back(static_cast<char>(event.text.unicode));
                }
            }
            // mouse press
            else if (event.type == sf::Event::MouseButtonPressed) {
                sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));

                // panel first
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
                                bodyPanel.visible = false;
                                selectedIndex = -1;
                                continue;
                            } else if (minRect.contains(mouseScreenF)) {
                                bodyPanel.minimized = !bodyPanel.minimized;
                                continue;
                            } else {
                                bodyPanel.dragging = true;
                                bodyPanel.dragOffset = mouseScreenF - bodyPanel.pos;
                                continue;
                            }
                        } else if (panelRect.contains(mouseScreenF)) {
                            continue; // click inside panel, ignore world
                        }
                    }
                }

                double cx = WIDTH * 0.5;
                double cy = HEIGHT * 0.5;
                Vec2 worldPos{
                    (mouseScreenF.x - cx) * metersPerPixel + viewCenterWorld.x,
                    (mouseScreenF.y - cy) * metersPerPixel + viewCenterWorld.y
                };

                if (event.mouseButton.button == sf::Mouse::Left) {
                    if (mode == InteractionMode::AddStill) {
                        Body b;
                        b.mass    = M_EARTH;
                        b.radius  = 6.0;
                        double vol = computeVolumeFromRadius(b.radius);
                        b.density  = b.mass / vol;
                        b.position = worldPos;
                        b.velocity = {0.0, 0.0};
                        b.color    = sf::Color::Cyan;
                        b.name     = "Body";
                        bodies.push_back(b);
                        computeGravity(bodies);
                    }
                    else if (mode == InteractionMode::AddMoving) {
                        isDraggingAddMoving = true;
                        dragStart = worldPos;
                    }
                    else if (mode == InteractionMode::AddOrbitingPlace) {
                        int parentIdx = findStrongestInfluenceParent(bodies, worldPos);
                        if (parentIdx < 0) {
                            std::cout << "No parent body found for orbiting placement.\n";
                        } else {
                            Body b;
                            b.mass    = M_EARTH;
                            b.radius  = 6.0;
                            double vol = computeVolumeFromRadius(b.radius);
                            b.density  = b.mass / vol;
                            b.position = worldPos;
                            b.color    = sf::Color::Green;
                            b.name     = "Orbiting body";
                            b.ghost    = true;

                            Body& parent = bodies[parentIdx];
                            Vec2 r = b.position - parent.position;
                            double dist = length(r);
                            if (dist < SOFTENING_METERS) dist = SOFTENING_METERS;

                            double vmag = std::sqrt(G_PHYS * parent.mass / dist);
                            Vec2 rhat{ r.x / dist, r.y / dist };
                            Vec2 perp{ -rhat.y, rhat.x };
                            Vec2 vOrb = perp * vmag;

                            b.velocity = parent.velocity + vOrb;

                            bodies.push_back(b);
                            orbitEditBodyIndex   = (int)bodies.size() - 1;
                            orbitEditParentIndex = parentIdx;
                            selectedIndex        = orbitEditBodyIndex;
                            bodyPanel.visible    = true;
                            bodyPanel.minimized  = false;
                            computeGravity(bodies); // ghosts ignored anyway
                            mode = InteractionMode::OrbitEdit;
                        }
                    }
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    double maxDistMeters = 25.0 * metersPerPixel;
                    int idx = findNearestBody(bodies, worldPos, maxDistMeters);
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
                            (mouseScreenF.x - cx) * metersPerPixel + viewCenterWorld.x,
                            (mouseScreenF.y - cy) * metersPerPixel + viewCenterWorld.y
                        };

                        Vec2 delta = worldEnd - dragStart;
                        Vec2 vel   = delta * 1e-5; // tweak factor

                        Body b;
                        b.mass    = M_EARTH;
                        b.radius  = 6.0;
                        double vol = computeVolumeFromRadius(b.radius);
                        b.density  = b.mass / vol;
                        b.position = dragStart;
                        b.velocity = vel;
                        b.color    = sf::Color::Green;
                        b.name     = "Body";
                        bodies.push_back(b);
                        computeGravity(bodies);
                    }
                }
            }
            // mouse move (for dragging panel)
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

        // AddMoving velocity arrow
        if (isDraggingAddMoving) {
            sf::Vector2f startScreen = worldToScreen(dragStart, WIDTH, HEIGHT);
            sf::Vector2i mousePix = sf::Mouse::getPosition(window);
            sf::Vector2f mouseScreenF = window.mapPixelToCoords(mousePix);
            sf::Vertex line[] = {
                sf::Vertex(startScreen, sf::Color::White),
                sf::Vertex(mouseScreenF, sf::Color::White)
            };
            window.draw(line, 2, sf::Lines);
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
            } else if (b.ghost) {
                circle.setOutlineThickness(1.5f);
                circle.setOutlineColor(sf::Color(0, 150, 255));
            }

            window.draw(circle);
        }

        // Orbit overlay in OrbitEdit
        if (fontLoaded &&
            mode == InteractionMode::OrbitEdit &&
            orbitEditBodyIndex >= 0 && orbitEditBodyIndex < (int)bodies.size() &&
            orbitEditParentIndex >= 0 && orbitEditParentIndex < (int)bodies.size()) {
            drawOrbitOverlay(window, bodies[orbitEditBodyIndex], bodies[orbitEditParentIndex],
                             WIDTH, HEIGHT, font);
        }

        // HUD
        if (fontLoaded) {
            sf::Text text;
            text.setFont(font);
            text.setCharacterSize(14);
            text.setFillColor(sf::Color::White);
            text.setPosition(10.0f, 10.0f);

            std::string modeStr;
            if (mode == InteractionMode::AddStill)       modeStr = "Add Still (1)";
            else if (mode == InteractionMode::AddMoving) modeStr = "Add Moving (2)";
            else if (mode == InteractionMode::AddOrbitingPlace || mode == InteractionMode::OrbitEdit)
                modeStr = "Add Orbiting (3)";

            std::string selStr = "None";
            if (selectedIndex >= 0 && selectedIndex < (int)bodies.size()) {
                selStr = bodies[selectedIndex].name;
            }

            std::string hud =
                "M1: add body   M2: select body\n"
                "1: add still   2: add moving   3: add orbiting (Enter=commit, Esc=cancel)\n"
                "Del: delete selected   Up/Down: change mass\n"
                "T: toggle trails  [ ]: trail length\n"
                "V: toggle volume lock   X: toggle fixed\n"
                "F: follow selected   F2: rename   Space: reset\n"
                "Z/X: zoom in/out   ,/. : slower/faster time\n" +
                std::string("Mode: ") + modeStr +
                "   Selected: " + selStr +
                (lockVolume ? "   [Volume LOCKED]" : "   [Volume UNLOCKED]");

            text.setString(hud);
            window.draw(text);

            // bottom-left: scale + time
            sf::Text scaleText;
            scaleText.setFont(font);
            scaleText.setCharacterSize(14);
            scaleText.setFillColor(sf::Color::White);
            scaleText.setPosition(10.0f, HEIGHT - 40.0f);
            scaleText.setString("Scale: " + formatScale());
            window.draw(scaleText);

            sf::Text timeText;
            timeText.setFont(font);
            timeText.setCharacterSize(14);
            timeText.setFillColor(sf::Color::White);
            timeText.setPosition(10.0f, HEIGHT - 22.0f);
            timeText.setString("Time:  " + formatTimeScaleStr());
            window.draw(timeText);
        }

        // Body panel
        if (fontLoaded && bodyPanel.visible && selectedIndex >= 0 && selectedIndex < (int)bodies.size()) {
            Body& b = bodies[selectedIndex];

            sf::RectangleShape panel;
            panel.setPosition(bodyPanel.pos);
            panel.setSize(bodyPanel.size);
            panel.setFillColor(sf::Color(20, 20, 35, 230));
            panel.setOutlineThickness(1.0f);
            panel.setOutlineColor(sf::Color::White);
            window.draw(panel);

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

                // Mass (scientific)
                line.setPosition(x, y);
                line.setString("Mass: " + formatSci(b.mass) + " kg");
                window.draw(line);
                y += 20.0f;

                // Density
                line.setPosition(x, y);
                line.setString("Density: " + formatSci(b.density));
                window.draw(line);
                y += 20.0f;

                // Radius (pixels)
                line.setPosition(x, y);
                line.setString("Radius: " + formatSci(b.radius) + " px");
                window.draw(line);
                y += 20.0f;

                // Flags
                std::string flags = "Fixed: ";
                flags += (b.fixed ? "Yes" : "No");
                flags += "   Volume lock: ";
                flags += (lockVolume ? "On" : "Off");
                line.setPosition(x, y);
                line.setString(flags);
                window.draw(line);
                y += 20.0f;
            }
        }

        window.display();
    }

    return 0;
}