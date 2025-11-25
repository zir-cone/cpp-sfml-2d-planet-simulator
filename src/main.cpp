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

double dot(const Vec2& a, const Vec2& b) {
    return a.x*b.x + a.y*b.y;
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
const double G_PHYS       = 6.67430e-11;          // m^3 kg^-1 s^-2
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

// formatting helpers
std::string formatSci(double v, int prec = 3) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*e", prec, v);
    return std::string(buf);
}

std::string formatDistance(double meters) {
    char buf[64];
    if (meters >= 1.0e12)      std::snprintf(buf, sizeof(buf), "%.3f Tm", meters / 1.0e12);
    else if (meters >= 1.0e9)  std::snprintf(buf, sizeof(buf), "%.3f Gm", meters / 1.0e9);
    else if (meters >= 1.0e6)  std::snprintf(buf, sizeof(buf), "%.3f Mm", meters / 1.0e6);
    else if (meters >= 1.0e3)  std::snprintf(buf, sizeof(buf), "%.3f km", meters / 1.0e3);
    else                       std::snprintf(buf, sizeof(buf), "%.0f m",  meters);
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

// --------------------
// orbit editing state
// --------------------
enum class OrbitHandle { None, Body, Pe, Ap };

struct OrbitEditState {
    bool active = false;
    int  bodyIndex = -1;
    int  parentIndex = -1;

    double a = 0.0;  // semi-major axis
    double e = 0.0;  // eccentricity (0..1)
    double nu = 0.0; // true anomaly (radians)

    Vec2 u{1.0, 0.0}; // unit vector from parent to periapsis
    Vec2 v{0.0, 1.0}; // perpendicular, CCW from u

    OrbitHandle dragHandle = OrbitHandle::None;
    bool dragging = false;
};

OrbitEditState orbitEdit;

// gui panel for sel body
struct BodyPanel {
    bool visible   = false;
    bool minimized = false;
    sf::Vector2f pos{20.0f, 80.0f};
    sf::Vector2f size{260.0f, 180.0f};
    bool dragging  = false;
    sf::Vector2f dragOffset{0.0f, 0.0f};
};

// reset confirmation popup
struct ResetPopup {
    bool visible = false;
    sf::FloatRect yesRect;
    sf::FloatRect noRect;
};

ResetPopup resetPopup;

// --------------------
// orbit helper functions
// --------------------
void normalize(Vec2& v) {
    double len = length(v);
    if (len > 0.0) {
        v.x /= len;
        v.y /= len;
    }
}

// update body.position & velocity from orbitEdit (around parent)
void updateBodyFromOrbit(Body& body, const Body& parent) {
    double a = orbitEdit.a;
    double e = orbitEdit.e;
    double nu = orbitEdit.nu;

    double cosnu = std::cos(nu);
    double sinnu = std::sin(nu);

    double denom = 1.0 + e * cosnu;
    if (denom <= 0.0) denom = 1e-6;
    double r = a * (1.0 - e * e) / denom;

    Vec2 rel = orbitEdit.u * (r * cosnu) + orbitEdit.v * (r * sinnu);
    body.position = { parent.position.x + rel.x, parent.position.y + rel.y };

    // vis-viva velocity magnitude
    double mu = G_PHYS * (body.mass + parent.mass);
    double vmag = std::sqrt(std::max(0.0, mu * (2.0 / r - 1.0 / a)));

    Vec2 rhat = rel;
    normalize(rhat);
    Vec2 that{ -rhat.y, rhat.x }; // CCW tangent

    body.velocity = { parent.velocity.x + that.x * vmag,
                      parent.velocity.y + that.y * vmag };
}

// compute Pe & Ap radii
double getRpe() {
    return orbitEdit.a * (1.0 - orbitEdit.e);
}
double getRap() {
    return orbitEdit.a * (1.0 + orbitEdit.e);
}

// init orbitEdit from current body & parent (start circular orbit)
void initOrbitEdit(int bodyIdx, int parentIdx, std::vector<Body>& bodies) {
    orbitEdit.active = true;
    orbitEdit.bodyIndex = bodyIdx;
    orbitEdit.parentIndex = parentIdx;
    orbitEdit.dragHandle = OrbitHandle::None;
    orbitEdit.dragging = false;

    Body& body = bodies[bodyIdx];
    Body& parent = bodies[parentIdx];

    Vec2 rel{ body.position.x - parent.position.x,
              body.position.y - parent.position.y };
    double r = length(rel);
    if (r < SOFTENING_METERS) r = SOFTENING_METERS;

    orbitEdit.a = r;
    orbitEdit.e = 0.0;
    orbitEdit.nu = 0.0;

    orbitEdit.u = rel;
    normalize(orbitEdit.u);
    orbitEdit.v = { -orbitEdit.u.y, orbitEdit.u.x };

    updateBodyFromOrbit(body, parent);
}

// orbit geometry for overlay & hit-testing
struct OrbitGeom {
    Vec2 parentPos;
    double a;
    double e;
    Vec2 u;
    Vec2 v;
    Vec2 peWorld;
    Vec2 apWorld;
};

bool computeOrbitGeom(const Body& parent, OrbitGeom& out) {
    if (!orbitEdit.active) return false;
    out.parentPos = parent.position;
    out.a = orbitEdit.a;
    out.e = orbitEdit.e;
    out.u = orbitEdit.u;
    out.v = orbitEdit.v;

    double r_pe = getRpe();
    double r_ap = getRap();
    if (r_pe <= 0.0 || r_ap <= 0.0) return false;

    out.peWorld = { parent.position.x + out.u.x * r_pe,
                    parent.position.y + out.u.y * r_pe };
    out.apWorld = { parent.position.x - out.u.x * r_ap,
                    parent.position.y - out.u.y * r_ap };
    return true;
}

// which overlay handle is near mouse?
OrbitHandle pickOrbitHandle(
    const Body& body,
    const Body& parent,
    const sf::Vector2f& mouseScreen,
    double width,
    double height
) {
    OrbitGeom geom;
    if (!computeOrbitGeom(parent, geom)) return OrbitHandle::None;

    sf::Vector2f peScreen = worldToScreen(geom.peWorld, width, height);
    sf::Vector2f apScreen = worldToScreen(geom.apWorld, width, height);
    sf::Vector2f bodyScreen = worldToScreen(body.position, width, height);

    auto dist2 = [](const sf::Vector2f& a, const sf::Vector2f& b) {
        float dx = a.x - b.x;
        float dy = a.y - b.y;
        return dx*dx + dy*dy;
    };

    const float r2 = 10.0f * 10.0f;

    if (dist2(mouseScreen, bodyScreen) <= r2) return OrbitHandle::Body;
    if (dist2(mouseScreen, peScreen)   <= r2) return OrbitHandle::Pe;
    if (dist2(mouseScreen, apScreen)   <= r2) return OrbitHandle::Ap;

    return OrbitHandle::None;
}

// apply dragging to orbit parameters
void applyOrbitDrag(
    Body& body,
    Body& parent,
    const Vec2& mouseWorld
) {
    if (!orbitEdit.active) return;

    Vec2 relMouse{
        mouseWorld.x - parent.position.x,
        mouseWorld.y - parent.position.y
    };
    double r_mouse = length(relMouse);
    if (r_mouse < 1.0e6) return;

    // Coordinates of relMouse in (u,v) basis
    double c = dot(relMouse, orbitEdit.u) / r_mouse;
    double s = dot(relMouse, orbitEdit.v) / r_mouse;
    double angleMouse = std::atan2(s, c);

    double r_pe_old = getRpe();
    double r_ap_old = getRap();

    if (orbitEdit.dragHandle == OrbitHandle::Body) {
        // move along ellipse (change ν, keep a,e)
        orbitEdit.nu = angleMouse;
    }
    else if (orbitEdit.dragHandle == OrbitHandle::Pe ||
             orbitEdit.dragHandle == OrbitHandle::Ap) {

        // new orientation: direction of mouse from parent
        orbitEdit.u = relMouse;
        normalize(orbitEdit.u);
        orbitEdit.v = { -orbitEdit.u.y, orbitEdit.u.x };

        double r_pe_new = r_pe_old;
        double r_ap_new = r_ap_old;

        if (orbitEdit.dragHandle == OrbitHandle::Pe) {
            // user intends to drag Pe distance
            r_pe_new = r_mouse;

            if (r_pe_new > r_ap_old) {
                // swap roles: mouse point is now Ap
                std::swap(r_pe_new, r_ap_old);
                orbitEdit.dragHandle = OrbitHandle::Ap;
                // orientation u points to new Pe, which is opposite current mouse direction:
                orbitEdit.u.x = -orbitEdit.u.x;
                orbitEdit.u.y = -orbitEdit.u.y;
                orbitEdit.v = { -orbitEdit.u.y, orbitEdit.u.x };
            }
        } else {
            // dragging Ap
            r_ap_new = r_mouse;

            if (r_ap_new < r_pe_old) {
                // swap roles: mouse point is now Pe
                std::swap(r_ap_new, r_pe_old);
                orbitEdit.dragHandle = OrbitHandle::Pe;
                // orientation u points toward new Pe (mouse direction)
                // (already set above)
            }
        }

        // recompute a,e from r_pe,r_ap
        r_pe_new = std::max(r_pe_new, 1.0e6);
        r_ap_new = std::max(r_ap_new, r_pe_new + 1.0e6);

        orbitEdit.a = 0.5 * (r_pe_new + r_ap_new);
        orbitEdit.e = (r_ap_new - r_pe_new) / (r_ap_new + r_pe_new);
        if (orbitEdit.e < 0.0) orbitEdit.e = 0.0;
        if (orbitEdit.e > 0.99) orbitEdit.e = 0.99;

        // keep ν as angleMouse (body stays roughly where it is angle-wise)
        orbitEdit.nu = angleMouse;
    }

    updateBodyFromOrbit(body, parent);
}

// --------------------
// drawing orbit overlay
// --------------------
void drawOrbitOverlay(
    sf::RenderWindow& window,
    const Body& body,
    const Body& parent,
    double width,
    double height,
    const sf::Font& font
) {
    if (!orbitEdit.active) return;

    OrbitGeom geom;
    if (!computeOrbitGeom(parent, geom)) return;

    // ellipse samples
    const int SEGMENTS = 256;
    std::vector<sf::Vertex> orbitLine;
    orbitLine.reserve(SEGMENTS + 1);

    for (int i = 0; i <= SEGMENTS; ++i) {
        double theta = (double)i / SEGMENTS * 2.0 * 3.141592653589793;
        double cosT = std::cos(theta);
        double sinT = std::sin(theta);
        double denom = 1.0 + geom.e * cosT;
        if (denom <= 0.0) denom = 1e-6;
        double r = geom.a * (1.0 - geom.e * geom.e) / denom;

        Vec2 rel = geom.u * (r * cosT) + geom.v * (r * sinT);
        Vec2 pWorld{ geom.parentPos.x + rel.x, geom.parentPos.y + rel.y };
        sf::Vector2f s = worldToScreen(pWorld, width, height);
        orbitLine.emplace_back(s, sf::Color(0, 150, 255));
    }
    window.draw(orbitLine.data(), orbitLine.size(), sf::LineStrip);

    // parent (central focus)
    sf::Vector2f parentScreen = worldToScreen(geom.parentPos, width, height);
    sf::CircleShape parentPt(4.0f);
    parentPt.setOrigin(4.0f, 4.0f);
    parentPt.setPosition(parentScreen);
    parentPt.setFillColor(sf::Color::Red);
    window.draw(parentPt);

    // Pe / Ap points
    sf::Vector2f peScreen = worldToScreen(geom.peWorld, width, height);
    sf::Vector2f apScreen = worldToScreen(geom.apWorld, width, height);

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

    // dashed lines parent->Pe/Ap
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

    drawDashed(parentScreen, peScreen);
    drawDashed(parentScreen, apScreen);

    if (font.getInfo().family.empty()) return;

    sf::Text text;
    text.setFont(font);
    text.setCharacterSize(12);

    // Pe / Ap labels
    text.setFillColor(sf::Color(255, 100, 200));
    text.setString("Pe");
    text.setPosition(peScreen.x + 6.0f, peScreen.y - 4.0f);
    window.draw(text);

    text.setString("Ap");
    text.setPosition(apScreen.x + 6.0f, apScreen.y - 4.0f);
    window.draw(text);

    // true anomaly ν and period
    Vec2 relBody{ body.position.x - geom.parentPos.x,
                  body.position.y - geom.parentPos.y };
    double r = length(relBody);
    if (r <= 0.0) r = 1.0;

    Vec2 rhat = relBody;
    normalize(rhat);
    double c = dot(rhat, geom.u);
    double s = dot(rhat, geom.v);
    double nuDeg = std::atan2(s, c) * 180.0 / 3.141592653589793;

    char buf[64];
    sf::Vector2f bodyScreen = worldToScreen(body.position, width, height);
    text.setFillColor(sf::Color::Green);
    std::snprintf(buf, sizeof(buf), "ν = %.2f°", nuDeg);
    text.setString(buf);
    text.setPosition(bodyScreen.x + 8.0f, bodyScreen.y - 20.0f);
    window.draw(text);

    double mu = G_PHYS * (body.mass + parent.mass);
    double T = 2.0 * 3.141592653589793 * std::sqrt(geom.a * geom.a * geom.a / mu);
    double T_years = T / YEAR_SECONDS;
    text.setFillColor(sf::Color(0, 180, 255));
    std::snprintf(buf, sizeof(buf), "P = %.3f yr", T_years);
    text.setString(buf);
    text.setPosition(parentScreen.x - 40.0f, parentScreen.y + 30.0f);
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

// --------------------
// reset system helper
// --------------------
void resetSystem(std::vector<Body>& bodies, Body& star) {
    bodies.clear();
    star.position = {0.0, 0.0};
    star.velocity = {0.0, 0.0};
    star.fixed    = true;
    bodies.push_back(star);
    computeGravity(bodies);
    for (auto& b : bodies) b.trail.clear();
    orbitEdit.active = false;
    orbitEdit.dragging = false;
    orbitEdit.bodyIndex = orbitEdit.parentIndex = -1;
    viewBodyIndex = -1;
    viewCenterWorld = star.position;
}

// --------------------
// main
// --------------------
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
                // Reset popup (Ctrl+R)
                if (event.key.code == sf::Keyboard::R && event.key.control) {
                    resetPopup.visible = true;
                }
                // cancel popup / orbit edit / rename via Esc
                else if (event.key.code == sf::Keyboard::Escape) {
                    if (resetPopup.visible) {
                        resetPopup.visible = false;
                    } else if (renaming) {
                        renaming = false;
                    } else if (mode == InteractionMode::OrbitEdit &&
                               orbitEdit.active &&
                               orbitEdit.bodyIndex >= 0 &&
                               orbitEdit.bodyIndex < (int)bodies.size()) {
                        // cancel ghost
                        bodies.erase(bodies.begin() + orbitEdit.bodyIndex);
                        selectedIndex = -1;
                        bodyPanel.visible = false;
                        orbitEdit.active = false;
                        orbitEdit.dragging = false;
                        computeGravity(bodies);
                        mode = InteractionMode::AddOrbitingPlace;
                    }
                }
                // orbit commit
                else if (event.key.code == sf::Keyboard::Enter ||
                         event.key.code == sf::Keyboard::Return) {
                    if (mode == InteractionMode::OrbitEdit &&
                        orbitEdit.active &&
                        orbitEdit.bodyIndex >= 0 &&
                        orbitEdit.bodyIndex < (int)bodies.size()) {
                        bodies[orbitEdit.bodyIndex].ghost = false;
                        computeGravity(bodies);
                        orbitEdit.active = false;
                        orbitEdit.dragging = false;
                        mode = InteractionMode::AddOrbitingPlace;
                    }
                }
                // other hotkeys only if popup not visible
                else if (!resetPopup.visible) {
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
                    }
                    // zoom I/O
                    else if (event.key.code == sf::Keyboard::I) { // zoom in
                        metersPerPixel *= 0.5;
                        if (metersPerPixel < 1.0e7) metersPerPixel = 1.0e7;
                    } else if (event.key.code == sf::Keyboard::O) { // zoom out
                        metersPerPixel *= 2.0;
                        if (metersPerPixel > 1.0e12) metersPerPixel = 1.0e12;
                    }
                    // time scale , / .
                    else if (event.key.code == sf::Keyboard::Comma) {
                        timeScale *= 0.5;
                        if (timeScale < 3600.0) timeScale = 3600.0;
                    } else if (event.key.code == sf::Keyboard::Period) {
                        timeScale *= 2.0;
                        if (timeScale > 1.0e7) timeScale = 1.0e7;
                    }
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

                // reset popup click
                if (resetPopup.visible && event.mouseButton.button == sf::Mouse::Left) {
                    if (resetPopup.yesRect.contains(mouseScreenF)) {
                        resetPopup.visible = false;
                        resetSystem(bodies, star);
                    } else if (resetPopup.noRect.contains(mouseScreenF)) {
                        resetPopup.visible = false;
                    }
                    continue;
                }

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

                // OrbitEdit handle picking
                if (mode == InteractionMode::OrbitEdit &&
                    event.mouseButton.button == sf::Mouse::Left &&
                    orbitEdit.active &&
                    orbitEdit.bodyIndex >= 0 && orbitEdit.bodyIndex < (int)bodies.size() &&
                    orbitEdit.parentIndex >= 0 && orbitEdit.parentIndex < (int)bodies.size()) {

                    OrbitHandle h = pickOrbitHandle(
                        bodies[orbitEdit.bodyIndex],
                        bodies[orbitEdit.parentIndex],
                        mouseScreenF,
                        WIDTH,
                        HEIGHT
                    );
                    if (h != OrbitHandle::None) {
                        orbitEdit.dragHandle = h;
                        orbitEdit.dragging = true;
                        continue;
                    }
                }

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
                            // initial circular guess
                            Vec2 r = { worldPos.x - parent.position.x,
                                       worldPos.y - parent.position.y };
                            double dist = length(r);
                            if (dist < SOFTENING_METERS) dist = SOFTENING_METERS;

                            double vmag = std::sqrt(G_PHYS * parent.mass / dist);
                            Vec2 rhat = r;
                            normalize(rhat);
                            Vec2 perp{ -rhat.y, rhat.x };
                            b.velocity = { parent.velocity.x + perp.x * vmag,
                                           parent.velocity.y + perp.y * vmag };

                            bodies.push_back(b);
                            int bodyIdx = (int)bodies.size() - 1;

                            initOrbitEdit(bodyIdx, parentIdx, bodies);
                            selectedIndex = bodyIdx;
                            bodyPanel.visible = true;
                            bodyPanel.minimized = false;
                            computeGravity(bodies); // ghost ignored anyway
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
                    if (orbitEdit.dragging) {
                        orbitEdit.dragging = false;
                        orbitEdit.dragHandle = OrbitHandle::None;
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
            // mouse move (panel + orbit drag)
            else if (event.type == sf::Event::MouseMoved) {
                if (bodyPanel.dragging) {
                    sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                    bodyPanel.pos = mouseScreenF - bodyPanel.dragOffset;
                }

                if (orbitEdit.dragging &&
                    orbitEdit.active &&
                    orbitEdit.bodyIndex >= 0 && orbitEdit.bodyIndex < (int)bodies.size() &&
                    orbitEdit.parentIndex >= 0 && orbitEdit.parentIndex < (int)bodies.size()) {

                    sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                    double cx = WIDTH * 0.5;
                    double cy = HEIGHT * 0.5;
                    Vec2 mouseWorld{
                        (mouseScreenF.x - cx) * metersPerPixel + viewCenterWorld.x,
                        (mouseScreenF.y - cy) * metersPerPixel + viewCenterWorld.y
                    };

                    applyOrbitDrag(
                        bodies[orbitEdit.bodyIndex],
                        bodies[orbitEdit.parentIndex],
                        mouseWorld
                    );
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
            orbitEdit.active &&
            orbitEdit.bodyIndex >= 0 && orbitEdit.bodyIndex < (int)bodies.size() &&
            orbitEdit.parentIndex >= 0 && orbitEdit.parentIndex < (int)bodies.size()) {
            drawOrbitOverlay(window, bodies[orbitEdit.bodyIndex], bodies[orbitEdit.parentIndex],
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
            else modeStr = "Add Orbiting (3)";

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
                "F: follow selected   F2: rename   Ctrl+R: reset (with confirm)\n"
                "I/O: zoom in/out   ,/. : slower/faster time\n" +
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

        // Reset confirmation popup
        if (resetPopup.visible && fontLoaded) {
            float pw = 300.0f;
            float ph = 140.0f;
            float px = (WIDTH - pw) * 0.5f;
            float py = (HEIGHT - ph) * 0.5f;

            sf::RectangleShape box(sf::Vector2f(pw, ph));
            box.setPosition(px, py);
            box.setFillColor(sf::Color(20, 20, 35, 240));
            box.setOutlineThickness(1.0f);
            box.setOutlineColor(sf::Color::White);
            window.draw(box);

            sf::Text msg;
            msg.setFont(font);
            msg.setCharacterSize(16);
            msg.setFillColor(sf::Color::White);
            msg.setString("Reset system to only the star?");
            msg.setPosition(px + 20.0f, py + 20.0f);
            window.draw(msg);

            float btnW = 80.0f;
            float btnH = 28.0f;
            float gap = 20.0f;

            float yesX = px + pw * 0.5f - gap - btnW;
            float noX  = px + pw * 0.5f + gap;
            float btnY = py + ph - btnH - 20.0f;

            sf::RectangleShape yesBtn(sf::Vector2f(btnW, btnH));
            yesBtn.setPosition(yesX, btnY);
            yesBtn.setFillColor(sf::Color(60, 120, 60));
            window.draw(yesBtn);

            sf::RectangleShape noBtn(sf::Vector2f(btnW, btnH));
            noBtn.setPosition(noX, btnY);
            noBtn.setFillColor(sf::Color(120, 60, 60));
            window.draw(noBtn);

            sf::Text yesTxt;
            yesTxt.setFont(font);
            yesTxt.setCharacterSize(14);
            yesTxt.setFillColor(sf::Color::White);
            yesTxt.setString("Yes");
            yesTxt.setPosition(yesX + 24.0f, btnY + 4.0f);
            window.draw(yesTxt);

            sf::Text noTxt;
            noTxt.setFont(font);
            noTxt.setCharacterSize(14);
            noTxt.setFillColor(sf::Color::White);
            noTxt.setString("No");
            noTxt.setPosition(noX + 28.0f, btnY + 4.0f);
            window.draw(noTxt);

            resetPopup.yesRect = sf::FloatRect(yesX, btnY, btnW, btnH);
            resetPopup.noRect  = sf::FloatRect(noX,  btnY, btnW, btnH);
        }

        window.display();
    }

    return 0;
}
