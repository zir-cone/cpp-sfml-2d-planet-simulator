#include <SFML/Graphics.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <cstdio>

// --------------------
// Basic 2d vector
// --------------------
struct Vec2 {
    double x, y;
    Vec2 operator+(const Vec2& o) const { return {x + o.x, y + o.y}; }
    Vec2 operator-(const Vec2& o) const { return {x - o.x, y - o.y}; }
    Vec2 operator*(double s)     const { return {x * s,   y * s};   }
    Vec2& operator+=(const Vec2& o) { x += o.x; y += o.y; return *this; }
    Vec2& operator*=(double s)      { x *= s;   y *= s;   return *this; }
};

double length(const Vec2& v) {
    return std::sqrt(v.x * v.x + v.y * v.y);
}
double dot(const Vec2& a, const Vec2& b) {
    return a.x * b.x + a.y * b.y;
}
void normalize(Vec2& v) {
    double len = length(v);
    if (len > 0.0) { v.x /= len; v.y /= len; }
}

// --------------------
// Formatting helpers
// --------------------
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
// Body types & composition
// --------------------
struct BodyTypeInfo {
    std::string name;
    std::string description;
};

enum class BodyCategory {
    Star,
    Rocky,
    Oceanic,
    GasGiant,
    IceGiant
};

struct Composition {
    double rock   = 0.0;
    double metal  = 0.0;
    double ice    = 0.0;
    double gas    = 0.0;
    double water  = 0.0;
};

// some default type infos
static BodyTypeInfo rockyType {
    "rocky terrestrial",
    "solid silicate/metal body with thin atmosphere"
};
static BodyTypeInfo oceanType {
    "oceanic terrestrial",
    "rocky body with large liquid-water oceans"
};
static BodyTypeInfo gasType {
    "gas giant",
    "massive hydrogen/helium-rich planet"
};
static BodyTypeInfo iceType {
    "ice giant",
    "planet rich in volatiles/ices"
};
static BodyTypeInfo starType {
    "main-sequence star",
    "self-luminous body powered by nuclear fusion"
};

// --------------------
// Bodies
// --------------------
struct Body {
    double mass;         // kg
    Vec2   position;     // m
    Vec2   velocity;     // m/s
    Vec2   acceleration; // m/s^2
    sf::Color color;

    double density;
    double radius;       // base visual radius at reference zoom (px)
    bool   fixed  = false;
    bool   ghost  = false;

    std::string  name;
    BodyTypeInfo typeInfo;
    BodyCategory category = BodyCategory::Rocky;
    Composition  comp;

    // thermal / radiative stuff
    double temperatureK   = 288.0;   // surface temp (planets) or effective temp (stars)
    double luminosityW    = 0.0;     // stars only
    double albedo         = 0.3;     // planets

    std::vector<Vec2> trail;
};

// --------------------
// Global constants
// --------------------
const double G_PHYS       = 6.67430e-11;              // m^3 kg^-1 s^-2
const double SIGMA_SB     = 5.670374419e-8;           // W m^-2 K^-4
const double YEAR_SECONDS = 365.25 * 24.0 * 3600.0;

const double M_SUN   = 1.98847e30;
const double M_EARTH = 5.97219e24;

// scale (mutable)
double metersPerPixel = 1.0e9;                        // 1 px = 1e9 m
double timeScale      = 86400.0 * 50.0;               // 50 days per real second

// softening
const double SOFTENING_METERS = 1.0e9;

bool trailsEnabled = true;
int  maxTrailPoints = 200;
bool lockVolume     = false;

// helpers for radius/density
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

// --------------------
// Camera
// --------------------
Vec2 viewCenterWorld{0.0, 0.0};
int  viewBodyIndex = -1;

sf::Vector2f worldToScreen(const Vec2& p, double width, double height) {
    double cx = width * 0.5;
    double cy = height * 0.5;

    double dx = (p.x - viewCenterWorld.x) / metersPerPixel;
    double dy = (p.y - viewCenterWorld.y) / metersPerPixel;

    return sf::Vector2f(
        static_cast<float>(cx + dx),
        static_cast<float>(cy + dy)
    );
}

// --------------------
// Gravity
// --------------------
void computeGravity(std::vector<Body>& bodies) {
    std::size_t n = bodies.size();
    for (auto& b : bodies) {
        b.acceleration = {0.0, 0.0};
    }

    for (std::size_t i = 0; i < n; ++i) {
        if (bodies[i].ghost) continue; // ghost: ignored
        for (std::size_t j = i + 1; j < n; ++j) {
            if (bodies[j].ghost) continue;

            Vec2 r = bodies[j].position - bodies[i].position;
            double dist2 = r.x * r.x + r.y * r.y;
            double softened2 = dist2 + SOFTENING_METERS * SOFTENING_METERS;
            double softened  = std::sqrt(softened2);
            if (softened < 1.0) softened = 1.0;
            double invDist3  = 1.0 / (softened2 * softened);

            Vec2 accel_i = r * (G_PHYS * bodies[j].mass * invDist3);
            Vec2 accel_j = r * (-G_PHYS * bodies[i].mass * invDist3);

            bodies[i].acceleration += accel_i;
            bodies[j].acceleration += accel_j;
        }
    }
}

// --------------------
// Thermal model
// --------------------
bool isStar(const Body& b) {
    return b.category == BodyCategory::Star;
}
void updateTemperatures(std::vector<Body>& bodies) {
    // collect stars
    std::vector<int> starIndices;
    for (int i = 0; i < (int)bodies.size(); ++i) {
        if (isStar(bodies[i]) && bodies[i].luminosityW > 0.0)
            starIndices.push_back(i);
    }
    if (starIndices.empty()) return;

    for (int i = 0; i < (int)bodies.size(); ++i) {
        Body& body = bodies[i];
        if (isStar(body)) continue;

        // find nearest star
        double bestDist2 = 1e99;
        const Body* bestStar = nullptr;
        for (int si : starIndices) {
            Body& s = bodies[si];
            Vec2 r = body.position - s.position;
            double d2 = r.x*r.x + r.y*r.y;
            if (d2 < bestDist2) {
                bestDist2 = d2;
                bestStar = &s;
            }
        }
        if (!bestStar) continue;

        double r = std::sqrt(bestDist2);
        if (r <= 0.0) r = 1.0;

        // simple radiative equilibrium
        double F = bestStar->luminosityW / (4.0 * 3.141592653589793 * r * r);
        double T = std::pow( (F * (1.0 - body.albedo)) / (4.0 * SIGMA_SB), 0.25 );
        if (std::isfinite(T)) body.temperatureK = T;
    }
}

// --------------------
// Step integrator
// --------------------
void stepSim(std::vector<Body>& bodies, double dtRealSeconds) {
    double dt = dtRealSeconds * timeScale;

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
    updateTemperatures(bodies);

    if (trailsEnabled) {
        for (auto& b : bodies) {
            if (!b.ghost) {
                b.trail.push_back(b.position);
                if ((int)b.trail.size() > maxTrailPoints)
                    b.trail.erase(b.trail.begin());
            }
        }
    } else {
        for (auto& b : bodies) b.trail.clear();
    }
}

// --------------------
// Interaction modes
// --------------------
enum class InteractionMode {
    AddStill,
    AddMoving,
    AddOrbitingPlace,
    OrbitEdit
};

int findNearestBody(const std::vector<Body>& bodies, const Vec2& pointMeters, double maxDistMeters) {
    int bestIndex = -1;
    double bestDist2 = maxDistMeters * maxDistMeters;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        double dx = bodies[i].position.x - pointMeters.x;
        double dy = bodies[i].position.y - pointMeters.y;
        double d2 = dx*dx + dy*dy;
        if (d2 < bestDist2) {
            bestDist2 = d2;
            bestIndex = (int)i;
        }
    }
    return bestIndex;
}
int findStrongestInfluenceParent(const std::vector<Body>& bodies, const Vec2& posMeters) {
    int bestIndex = -1;
    double bestAccel = 0.0;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        if (bodies[i].ghost) continue;
        Vec2 r = bodies[i].position - posMeters;
        double dist2 = r.x*r.x + r.y*r.y;
        if (dist2 <= 0.0) continue;
        double accel = G_PHYS * bodies[i].mass / dist2;
        if (accel > bestAccel) {
            bestAccel = accel;
            bestIndex = (int)i;
        }
    }
    return bestIndex;
}

// --------------------
// Panels (body editor, controls)
// --------------------
struct Panel {
    bool visible   = false;
    bool minimized = false;
    sf::Vector2f pos{20.f, 80.f};
    sf::Vector2f size{260.f, 180.f};
    bool dragging  = false;
    sf::Vector2f dragOffset{0.f, 0.f};
};
bool pointInRect(const sf::Vector2f& p, const sf::FloatRect& r) {
    return r.contains(p);
}

// --------------------
// Orbit edit
// --------------------
enum class OrbitHandle { None, Body, Pe, Ap };

struct OrbitEditState {
    bool active = false;
    int  bodyIndex   = -1;
    int  parentIndex = -1;

    double a  = 0.0;
    double e  = 0.0;
    double nu = 0.0;

    Vec2 u{1.0, 0.0};
    Vec2 v{0.0, 1.0};

    OrbitHandle dragHandle = OrbitHandle::None;
    bool dragging = false;
};

OrbitEditState orbitEdit;

double getRpe() { return orbitEdit.a * (1.0 - orbitEdit.e); }
double getRap() { return orbitEdit.a * (1.0 + orbitEdit.e); }

void updateBodyFromOrbit(Body& body, const Body& parent) {
    double a  = orbitEdit.a;
    double e  = orbitEdit.e;
    double nu = orbitEdit.nu;

    double cosnu = std::cos(nu);
    double sinnu = std::sin(nu);
    double denom = 1.0 + e * cosnu;
    if (denom <= 0.0) denom = 1e-6;
    double r = a * (1.0 - e*e) / denom;

    Vec2 rel = orbitEdit.u * (r * cosnu) + orbitEdit.v * (r * sinnu);
    body.position = { parent.position.x + rel.x, parent.position.y + rel.y };

    double mu = G_PHYS * (body.mass + parent.mass);
    double vmag = std::sqrt(std::max(0.0, mu * (2.0 / r - 1.0 / a)));
    Vec2 rhat = rel; normalize(rhat);
    Vec2 that{ -rhat.y, rhat.x };

    body.velocity = { parent.velocity.x + that.x * vmag,
                      parent.velocity.y + that.y * vmag };
}

void initOrbitEdit(int bodyIdx, int parentIdx, std::vector<Body>& bodies) {
    orbitEdit.active      = true;
    orbitEdit.bodyIndex   = bodyIdx;
    orbitEdit.parentIndex = parentIdx;
    orbitEdit.dragging    = false;
    orbitEdit.dragHandle  = OrbitHandle::None;

    Body& body   = bodies[bodyIdx];
    Body& parent = bodies[parentIdx];

    Vec2 rel{ body.position.x - parent.position.x,
              body.position.y - parent.position.y };
    double r = length(rel);
    if (r < SOFTENING_METERS) r = SOFTENING_METERS;

    orbitEdit.a  = r;
    orbitEdit.e  = 0.0;
    orbitEdit.nu = 0.0;

    orbitEdit.u = rel; normalize(orbitEdit.u);
    orbitEdit.v = { -orbitEdit.u.y, orbitEdit.u.x };

    updateBodyFromOrbit(body, parent);
}

struct OrbitGeom {
    Vec2 parentPos;
    double a;
    double e;
    Vec2 u, v;
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

OrbitHandle pickOrbitHandle(
    const Body& body,
    const Body& parent,
    const sf::Vector2f& mouseScreen,
    double width,
    double height
) {
    OrbitGeom geom;
    if (!computeOrbitGeom(parent, geom)) return OrbitHandle::None;

    sf::Vector2f peScreen   = worldToScreen(geom.peWorld, width, height);
    sf::Vector2f apScreen   = worldToScreen(geom.apWorld, width, height);
    sf::Vector2f bodyScreen = worldToScreen(body.position, width, height);

    auto d2 = [](const sf::Vector2f& a, const sf::Vector2f& b) {
        float dx = a.x - b.x, dy = a.y - b.y;
        return dx*dx + dy*dy;
    };

    const float r2 = 10.f * 10.f;
    if (d2(mouseScreen, bodyScreen) <= r2) return OrbitHandle::Body;
    if (d2(mouseScreen, peScreen)   <= r2) return OrbitHandle::Pe;
    if (d2(mouseScreen, apScreen)   <= r2) return OrbitHandle::Ap;
    return OrbitHandle::None;
}

void applyOrbitDrag(Body& body, Body& parent, const Vec2& mouseWorld) {
    if (!orbitEdit.active) return;

    Vec2 relMouse{ mouseWorld.x - parent.position.x,
                   mouseWorld.y - parent.position.y };
    double r_mouse = length(relMouse);
    if (r_mouse < 1.0e6) return;

    double c = dot(relMouse, orbitEdit.u) / r_mouse;
    double s = dot(relMouse, orbitEdit.v) / r_mouse;
    double angleMouse = std::atan2(s, c);

    double r_pe_old = getRpe();
    double r_ap_old = getRap();

    if (orbitEdit.dragHandle == OrbitHandle::Body) {
        orbitEdit.nu = angleMouse;
    } else if (orbitEdit.dragHandle == OrbitHandle::Pe ||
               orbitEdit.dragHandle == OrbitHandle::Ap) {
        orbitEdit.u = relMouse; normalize(orbitEdit.u);
        orbitEdit.v = { -orbitEdit.u.y, orbitEdit.u.x };

        double r_pe_new = r_pe_old;
        double r_ap_new = r_ap_old;

        if (orbitEdit.dragHandle == OrbitHandle::Pe) {
            r_pe_new = r_mouse;
            if (r_pe_new > r_ap_old) {
                std::swap(r_pe_new, r_ap_old);
                orbitEdit.dragHandle = OrbitHandle::Ap;
                orbitEdit.u.x = -orbitEdit.u.x;
                orbitEdit.u.y = -orbitEdit.u.y;
                orbitEdit.v   = { -orbitEdit.u.y, orbitEdit.u.x };
            }
        } else {
            r_ap_new = r_mouse;
            if (r_ap_new < r_pe_old) {
                std::swap(r_ap_new, r_pe_old);
                orbitEdit.dragHandle = OrbitHandle::Pe;
            }
        }

        r_pe_new = std::max(r_pe_new, 1.0e6);
        r_ap_new = std::max(r_ap_new, r_pe_new + 1.0e6);

        orbitEdit.a = 0.5 * (r_pe_new + r_ap_new);
        orbitEdit.e = (r_ap_new - r_pe_new) / (r_ap_new + r_pe_new);
        if (orbitEdit.e < 0.0)   orbitEdit.e = 0.0;
        if (orbitEdit.e > 0.99)  orbitEdit.e = 0.99;

        orbitEdit.nu = angleMouse;
    }

    updateBodyFromOrbit(body, parent);
}

// --------------------
// Colour functions
// --------------------
sf::Color colorFromStarTemp(double T) {
    // very crude Kelvin -> RGB approximation
    T = std::clamp(T, 2000.0, 40000.0);
    double t = T / 100.0;

    double r, g, b;
    // based on Tanner Helland's approximation
    if (t <= 66.0) {
        r = 255.0;
        g = 99.4708025861 * std::log(t) - 161.1195681661;
        if (t <= 19.0)
            b = 0.0;
        else {
            b = 138.5177312231 * std::log(t - 10.0) - 305.0447927307;
        }
    } else {
        r = 329.698727446 * std::pow(t - 60.0, -0.1332047592);
        g = 288.1221695283 * std::pow(t - 60.0, -0.0755148492);
        b = 255.0;
    }
    auto clamp255 = [](double x) {
        if (x < 0.0) x = 0.0;
        if (x > 255.0) x = 255.0;
        return (sf::Uint8)x;
    };
    return sf::Color(clamp255(r), clamp255(g), clamp255(b));
}

sf::Color colorFromPlanet(const Body& b) {
    // base hues by category / composition
    sf::Color base;

    if (b.category == BodyCategory::Rocky) {
        base = sf::Color(150, 120, 90);
    } else if (b.category == BodyCategory::Oceanic) {
        base = sf::Color(60, 110, 170);
    } else if (b.category == BodyCategory::GasGiant) {
        base = sf::Color(200, 170, 120);
    } else if (b.category == BodyCategory::IceGiant) {
        base = sf::Color(100, 150, 220);
    } else {
        base = sf::Color(150, 150, 150);
    }

    // modulate brightness with temperature (warmer -> slightly brighter)
    double T = b.temperatureK;
    double factor = std::clamp((T - 100.0) / 400.0, 0.6, 1.4);
    auto cl = [&](int c) {
        double v = c * factor;
        if (v < 0.0) v = 0.0;
        if (v > 255.0) v = 255.0;
        return (sf::Uint8)v;
    };
    return sf::Color(cl(base.r), cl(base.g), cl(base.b));
}

sf::Color bodyColor(const Body& b) {
    if (isStar(b)) {
        return colorFromStarTemp(b.temperatureK);
    } else {
        return colorFromPlanet(b);
    }
}

// --------------------
// Reset popup
// --------------------
struct ResetPopup {
    bool visible = false;
    sf::FloatRect yesRect;
    sf::FloatRect noRect;
};
ResetPopup resetPopup;

// scale/time formatting
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
// Reset system helper
// --------------------
void resetSystem(std::vector<Body>& bodies, Body& star) {
    bodies.clear();
    star.position = {0.0, 0.0};
    star.velocity = {0.0, 0.0};
    star.fixed    = true;
    star.ghost    = false;
    bodies.push_back(star);
    computeGravity(bodies);
    updateTemperatures(bodies);
    for (auto& b : bodies) b.trail.clear();
    orbitEdit = {};
    viewBodyIndex = -1;
    viewCenterWorld = star.position;
}

// --------------------
// Draw dashed line helper
// --------------------
void drawDashedLine(sf::RenderWindow& window,
                    const sf::Vector2f& A,
                    const sf::Vector2f& B,
                    sf::Color col = sf::Color::White) {
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
            sf::Vertex(s, col),
            sf::Vertex(e, col)
        };
        window.draw(line, 2, sf::Lines);
        pos += dashLen + gapLen;
    }
}

// --------------------
// Orbit overlay drawing
// --------------------
void drawOrbitOverlay(
    sf::RenderWindow& window,
    const Body& body,
    const Body& parent,
    double WIDTH,
    double HEIGHT,
    const sf::Font& font
) {
    if (!orbitEdit.active) return;

    OrbitGeom geom;
    if (!computeOrbitGeom(parent, geom)) return;

    // ellipse line
    const int SEG = 256;
    std::vector<sf::Vertex> orbitLine;
    orbitLine.reserve(SEG + 1);

    for (int i = 0; i <= SEG; ++i) {
        double th = (double)i / SEG * 2.0 * 3.141592653589793;
        double cosT = std::cos(th);
        double sinT = std::sin(th);
        double denom = 1.0 + geom.e * cosT;
        if (denom <= 0.0) denom = 1e-6;
        double r = geom.a * (1.0 - geom.e*geom.e) / denom;
        Vec2 rel = geom.u * (r * cosT) + geom.v * (r * sinT);
        Vec2 pWorld{ geom.parentPos.x + rel.x, geom.parentPos.y + rel.y };
        orbitLine.emplace_back(worldToScreen(pWorld, WIDTH, HEIGHT),
                               sf::Color(0, 150, 255));
    }
    window.draw(orbitLine.data(), orbitLine.size(), sf::LineStrip);

    sf::Vector2f parentScreen = worldToScreen(geom.parentPos, WIDTH, HEIGHT);
    sf::CircleShape parentPt(4.f);
    parentPt.setOrigin(4.f, 4.f);
    parentPt.setPosition(parentScreen);
    parentPt.setFillColor(sf::Color::Red);
    window.draw(parentPt);

    sf::Vector2f peScreen = worldToScreen(geom.peWorld, WIDTH, HEIGHT);
    sf::Vector2f apScreen = worldToScreen(geom.apWorld, WIDTH, HEIGHT);

    sf::CircleShape pePt(4.f);
    pePt.setOrigin(4.f, 4.f);
    pePt.setPosition(peScreen);
    pePt.setFillColor(sf::Color(255, 100, 200));
    window.draw(pePt);

    sf::CircleShape apPt(4.f);
    apPt.setOrigin(4.f, 4.f);
    apPt.setPosition(apScreen);
    apPt.setFillColor(sf::Color(255, 100, 200));
    window.draw(apPt);

    drawDashedLine(window, parentScreen, peScreen);
    drawDashedLine(window, parentScreen, apScreen);

    if (font.getInfo().family.empty()) return;

    sf::Text text;
    text.setFont(font);
    text.setCharacterSize(12);

    // Pe/Ap labels
    text.setFillColor(sf::Color(255, 100, 200));
    text.setString("Pe");
    text.setPosition(peScreen.x + 6.f, peScreen.y - 4.f);
    window.draw(text);

    text.setString("Ap");
    text.setPosition(apScreen.x + 6.f, apScreen.y - 4.f);
    window.draw(text);

    // distances at Pe / Ap (from center)
    double rPe = getRpe();
    double rAp = getRap();
    sf::Vector2f midPe = (parentScreen + peScreen) * 0.5f;
    sf::Vector2f midAp = (parentScreen + apScreen) * 0.5f;

    text.setFillColor(sf::Color::White);
    text.setPosition(midPe.x + 4.f, midPe.y - 12.f);
    text.setString("r(Pe) = " + formatDistance(rPe));
    window.draw(text);

    text.setPosition(midAp.x + 4.f, midAp.y - 12.f);
    text.setString("r(Ap) = " + formatDistance(rAp));
    window.draw(text);

    // true anomaly & period
    Vec2 relBody{ body.position.x - geom.parentPos.x,
                  body.position.y - geom.parentPos.y };
    double r = length(relBody); if (r <= 0.0) r = 1.0;
    Vec2 rhat = relBody; normalize(rhat);
    double c = dot(rhat, geom.u);
    double s = dot(rhat, geom.v);
    double nuDeg = std::atan2(s, c) * 180.0 / 3.141592653589793;

    char buf[64];
    sf::Vector2f bodyScreen = worldToScreen(body.position, WIDTH, HEIGHT);
    text.setFillColor(sf::Color::Green);
    std::snprintf(buf, sizeof(buf), "ν = %.2f°", nuDeg);
    text.setString(buf);
    text.setPosition(bodyScreen.x + 8.f, bodyScreen.y - 20.f);
    window.draw(text);

    double mu = G_PHYS * (body.mass + parent.mass);
    double T = 2.0 * 3.141592653589793 * std::sqrt(geom.a*geom.a*geom.a / mu);
    double T_years = T / YEAR_SECONDS;
    text.setFillColor(sf::Color(0, 180, 255));
    std::snprintf(buf, sizeof(buf), "P = %.3f yr", T_years);
    text.setString(buf);
    text.setPosition(parentScreen.x - 40.f, parentScreen.y + 30.f);
    window.draw(text);
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

    // star
    Body star;
    star.mass        = M_SUN;
    star.position    = {0.0, 0.0};
    star.velocity    = {0.0, 0.0};
    star.fixed       = true;
    star.density     = 1.0;
    star.radius      = 18.0;
    star.name        = "Sol";
    star.category    = BodyCategory::Star;
    star.typeInfo    = starType;
    star.temperatureK = 5800.0;
    star.luminosityW  = 3.828e26;
    star.comp.gas     = 0.9;
    star.comp.metal   = 0.1;
    star.albedo       = 0.0;
    star.color        = bodyColor(star);

    bodies.push_back(star);
    computeGravity(bodies);
    updateTemperatures(bodies);

    InteractionMode mode = InteractionMode::AddStill;
    int selectedIndex = -1;

    bool isDraggingAddMoving = false;
    Vec2 dragStart{0.0, 0.0};

    sf::Clock clock;
    sf::Font font;
    bool fontLoaded = font.loadFromFile("resources/arial.ttf");

    Panel bodyPanel;
    bodyPanel.pos  = {20.f, 80.f};
    bodyPanel.size = {280.f, 210.f};

    Panel controlsPanel;
    controlsPanel.pos  = {WIDTH - 380.f, 10.f};
    controlsPanel.size = {360.f, 210.f};
    controlsPanel.visible = true;

    bool renaming = false;
    std::string renameBuffer;

    viewCenterWorld = star.position;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            // ------------- events -------------
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed) {
                // Ctrl+R reset popup
                if (event.key.code == sf::Keyboard::R && event.key.control) {
                    resetPopup.visible = true;
                }
                // toggle controls panel with '?'
                else if (event.key.code == sf::Keyboard::Slash && event.key.shift) {
                    controlsPanel.visible = !controlsPanel.visible;
                }
                // Esc: close dialogs / rename / orbit edit
                else if (event.key.code == sf::Keyboard::Escape) {
                    if (resetPopup.visible) {
                        resetPopup.visible = false;
                    } else if (renaming) {
                        renaming = false;
                    } else if (mode == InteractionMode::OrbitEdit &&
                               orbitEdit.active &&
                               orbitEdit.bodyIndex >= 0 &&
                               orbitEdit.bodyIndex < (int)bodies.size()) {
                        bodies.erase(bodies.begin() + orbitEdit.bodyIndex);
                        selectedIndex = -1;
                        bodyPanel.visible = false;
                        orbitEdit = {};
                        computeGravity(bodies);
                        mode = InteractionMode::AddOrbitingPlace;
                    }
                }
                // Enter to commit orbit
                else if ((event.key.code == sf::Keyboard::Enter ||
                          event.key.code == sf::Keyboard::Return) &&
                         mode == InteractionMode::OrbitEdit &&
                         orbitEdit.active &&
                         orbitEdit.bodyIndex >= 0 &&
                         orbitEdit.bodyIndex < (int)bodies.size()) {
                    bodies[orbitEdit.bodyIndex].ghost = false;
                    computeGravity(bodies);
                    updateTemperatures(bodies);
                    orbitEdit = {};
                    mode = InteractionMode::AddOrbitingPlace;
                }
                else if (!resetPopup.visible) {
                    // mode keys
                    if (event.key.code == sf::Keyboard::Num1)
                        mode = InteractionMode::AddStill;
                    else if (event.key.code == sf::Keyboard::Num2)
                        mode = InteractionMode::AddMoving;
                    else if (event.key.code == sf::Keyboard::Num3)
                        mode = InteractionMode::AddOrbitingPlace;
                    else if (event.key.code == sf::Keyboard::Delete) {
                        if (selectedIndex >= 0 && selectedIndex < (int)bodies.size()) {
                            bodies.erase(bodies.begin() + selectedIndex);
                            selectedIndex = -1;
                            bodyPanel.visible = false;
                            computeGravity(bodies);
                        }
                    }
                    else if (event.key.code == sf::Keyboard::Up) {
                        if (selectedIndex >= 0) {
                            bodies[selectedIndex].mass *= 1.2;
                            updateBodyAfterMassChange(bodies[selectedIndex]);
                        }
                    }
                    else if (event.key.code == sf::Keyboard::Down) {
                        if (selectedIndex >= 0) {
                            bodies[selectedIndex].mass *= 0.8;
                            updateBodyAfterMassChange(bodies[selectedIndex]);
                        }
                    }
                    else if (event.key.code == sf::Keyboard::T) {
                        trailsEnabled = !trailsEnabled;
                    }
                    else if (event.key.code == sf::Keyboard::LBracket) {
                        if (maxTrailPoints > 10) maxTrailPoints -= 10;
                    }
                    else if (event.key.code == sf::Keyboard::RBracket) {
                        maxTrailPoints += 10;
                    }
                    else if (event.key.code == sf::Keyboard::V) {
                        lockVolume = !lockVolume;
                    }
                    else if (event.key.code == sf::Keyboard::F) {
                        if (selectedIndex >= 0) {
                            if (viewBodyIndex == selectedIndex) viewBodyIndex = -1;
                            else viewBodyIndex = selectedIndex;
                        }
                    }
                    else if (event.key.code == sf::Keyboard::X) {
                        if (selectedIndex >= 0) {
                            bodies[selectedIndex].fixed = !bodies[selectedIndex].fixed;
                        }
                    }
                    else if (event.key.code == sf::Keyboard::F2) {
                        if (selectedIndex >= 0) {
                            renaming = true;
                            renameBuffer = bodies[selectedIndex].name;
                        }
                    }
                    // zoom + time
                    else if (event.key.code == sf::Keyboard::I) {
                        metersPerPixel *= 0.5;
                        if (metersPerPixel < 1.0e7) metersPerPixel = 1.0e7;
                    }
                    else if (event.key.code == sf::Keyboard::O) {
                        metersPerPixel *= 2.0;
                        if (metersPerPixel > 1.0e12) metersPerPixel = 1.0e12;
                    }
                    else if (event.key.code == sf::Keyboard::Comma) {
                        timeScale *= 0.5;
                        if (timeScale < 3600.0) timeScale = 3600.0;
                    }
                    else if (event.key.code == sf::Keyboard::Period) {
                        timeScale *= 2.0;
                        if (timeScale > 1.0e8) timeScale = 1.0e8;
                    }
                }
            }
            // rename text
            else if (event.type == sf::Event::TextEntered && renaming) {
                if (event.text.unicode == '\r' || event.text.unicode == '\n') {
                    if (selectedIndex >= 0) bodies[selectedIndex].name = renameBuffer;
                    renaming = false;
                } else if (event.text.unicode == 8) {
                    if (!renameBuffer.empty()) renameBuffer.pop_back();
                } else if (event.text.unicode >= 32 && event.text.unicode < 127) {
                    renameBuffer.push_back(static_cast<char>(event.text.unicode));
                }
            }
            // mouse pressed
            else if (event.type == sf::Event::MouseButtonPressed) {
                sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));

                // reset popup
                if (resetPopup.visible && event.mouseButton.button == sf::Mouse::Left) {
                    if (resetPopup.yesRect.contains(mouseScreenF)) {
                        resetPopup.visible = false;
                        resetSystem(bodies, star);
                    } else if (resetPopup.noRect.contains(mouseScreenF)) {
                        resetPopup.visible = false;
                    }
                    continue;
                }

                // panel handling (body)
                auto handlePanelClick = [&](Panel& p, bool allowClose) -> bool {
                    if (!p.visible) return false;
                    sf::FloatRect panelRect(p.pos, p.size);
                    sf::FloatRect titleRect(p.pos.x, p.pos.y, p.size.x, 24.f);

                    if (event.mouseButton.button == sf::Mouse::Left) {
                        if (titleRect.contains(mouseScreenF)) {
                            float btn = 16.f;
                            sf::FloatRect closeRect(
                                p.pos.x + p.size.x - btn - 4.f,
                                p.pos.y + 4.f, btn, btn);
                            sf::FloatRect minRect(
                                p.pos.x + p.size.x - 2.f*(btn+4.f),
                                p.pos.y + 4.f, btn, btn);
                            if (allowClose && closeRect.contains(mouseScreenF)) {
                                p.visible = false;
                                return true;
                            } else if (minRect.contains(mouseScreenF)) {
                                p.minimized = !p.minimized;
                                return true;
                            } else {
                                p.dragging = true;
                                p.dragOffset = mouseScreenF - p.pos;
                                return true;
                            }
                        } else if (panelRect.contains(mouseScreenF)) {
                            return true; // click absorbed
                        }
                    }
                    return false;
                };

                if (handlePanelClick(bodyPanel, true)) continue;
                if (handlePanelClick(controlsPanel, false)) continue;

                double cx = WIDTH * 0.5;
                double cy = HEIGHT * 0.5;
                Vec2 worldPos{
                    (mouseScreenF.x - cx) * metersPerPixel + viewCenterWorld.x,
                    (mouseScreenF.y - cy) * metersPerPixel + viewCenterWorld.y
                };

                // orbit handle pick
                if (mode == InteractionMode::OrbitEdit &&
                    event.mouseButton.button == sf::Mouse::Left &&
                    orbitEdit.active &&
                    orbitEdit.bodyIndex >= 0 &&
                    orbitEdit.bodyIndex < (int)bodies.size() &&
                    orbitEdit.parentIndex >= 0 &&
                    orbitEdit.parentIndex < (int)bodies.size()) {
                    OrbitHandle h = pickOrbitHandle(
                        bodies[orbitEdit.bodyIndex],
                        bodies[orbitEdit.parentIndex],
                        mouseScreenF, WIDTH, HEIGHT);
                    if (h != OrbitHandle::None) {
                        orbitEdit.dragHandle = h;
                        orbitEdit.dragging   = true;
                        continue;
                    }
                }

                if (event.mouseButton.button == sf::Mouse::Left) {
                    if (mode == InteractionMode::AddStill) {
                        Body b;
                        b.mass     = M_EARTH;
                        b.radius   = 6.0;
                        double vol = computeVolumeFromRadius(b.radius);
                        b.density  = b.mass / vol;
                        b.position = worldPos;
                        b.velocity = {0.0, 0.0};
                        b.category = BodyCategory::Rocky;
                        b.typeInfo = rockyType;
                        b.comp.rock = 0.7; b.comp.metal = 0.3;
                        b.albedo   = 0.3;
                        b.name     = "Rocky body";
                        b.color    = bodyColor(b);
                        bodies.push_back(b);
                        computeGravity(bodies);
                        updateTemperatures(bodies);
                    }
                    else if (mode == InteractionMode::AddMoving) {
                        isDraggingAddMoving = true;
                        dragStart = worldPos;
                    }
                    else if (mode == InteractionMode::AddOrbitingPlace) {
                        int parentIdx = findStrongestInfluenceParent(bodies, worldPos);
                        if (parentIdx >= 0) {
                            Body b;
                            b.mass     = M_EARTH;
                            b.radius   = 6.0;
                            double vol = computeVolumeFromRadius(b.radius);
                            b.density  = b.mass / vol;
                            b.position = worldPos;
                            b.category = BodyCategory::Rocky;
                            b.typeInfo = rockyType;
                            b.comp.rock = 0.7; b.comp.metal = 0.3;
                            b.albedo   = 0.3;
                            b.name     = "Orbiting body";
                            b.ghost    = true;

                            Body& parent = bodies[parentIdx];
                            Vec2 r = { worldPos.x - parent.position.x,
                                       worldPos.y - parent.position.y };
                            double dist = length(r);
                            if (dist < SOFTENING_METERS) dist = SOFTENING_METERS;
                            double vmag = std::sqrt(G_PHYS * parent.mass / dist);
                            Vec2 rhat = r; normalize(rhat);
                            Vec2 perp{ -rhat.y, rhat.x };
                            b.velocity = { parent.velocity.x + perp.x * vmag,
                                           parent.velocity.y + perp.y * vmag };
                            b.color = bodyColor(b);

                            bodies.push_back(b);
                            int bodyIdx = (int)bodies.size() - 1;
                            initOrbitEdit(bodyIdx, parentIdx, bodies);
                            selectedIndex = bodyIdx;
                            bodyPanel.visible = true;
                            bodyPanel.minimized = false;
                            mode = InteractionMode::OrbitEdit;
                        }
                    }
                } else if (event.mouseButton.button == sf::Mouse::Right) {
                    double maxDistMeters = 25.0 * metersPerPixel;
                    int idx = findNearestBody(bodies, worldPos, maxDistMeters);
                    selectedIndex = idx;
                    if (idx >= 0) {
                        bodyPanel.visible   = true;
                        bodyPanel.minimized = false;
                    }
                }
            }
            // mouse released
            else if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    bodyPanel.dragging = false;
                    controlsPanel.dragging = false;
                    if (orbitEdit.dragging) {
                        orbitEdit.dragging   = false;
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
                        Vec2 vel   = delta * 1e-5;

                        Body b;
                        b.mass     = M_EARTH;
                        b.radius   = 6.0;
                        double vol = computeVolumeFromRadius(b.radius);
                        b.density  = b.mass / vol;
                        b.position = dragStart;
                        b.velocity = vel;
                        b.category = BodyCategory::Rocky;
                        b.typeInfo = rockyType;
                        b.comp.rock = 0.7; b.comp.metal = 0.3;
                        b.albedo   = 0.3;
                        b.name     = "Moving body";
                        b.color    = bodyColor(b);
                        bodies.push_back(b);
                        computeGravity(bodies);
                        updateTemperatures(bodies);
                    }
                }
            }
            // mouse moved
            else if (event.type == sf::Event::MouseMoved) {
                sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                if (bodyPanel.dragging) {
                    bodyPanel.pos = mouseScreenF - bodyPanel.dragOffset;
                }
                if (controlsPanel.dragging) {
                    controlsPanel.pos = mouseScreenF - controlsPanel.dragOffset;
                }
                if (orbitEdit.dragging &&
                    orbitEdit.active &&
                    orbitEdit.bodyIndex >= 0 &&
                    orbitEdit.bodyIndex < (int)bodies.size() &&
                    orbitEdit.parentIndex >= 0 &&
                    orbitEdit.parentIndex < (int)bodies.size()) {

                    double cx = WIDTH * 0.5;
                    double cy = HEIGHT * 0.5;
                    Vec2 mouseWorld{
                        (mouseScreenF.x - cx) * metersPerPixel + viewCenterWorld.x,
                        (mouseScreenF.y - cy) * metersPerPixel + viewCenterWorld.y
                    };
                    applyOrbitDrag(bodies[orbitEdit.bodyIndex], bodies[orbitEdit.parentIndex], mouseWorld);
                }
            }
        } // end poll

        // -------- update physics --------
        float dtSeconds = clock.restart().asSeconds();
        if (dtSeconds > 0.05f) dtSeconds = 0.05f;
        stepSim(bodies, dtSeconds);

        if (viewBodyIndex >= 0 && viewBodyIndex < (int)bodies.size()) {
            viewCenterWorld = bodies[viewBodyIndex].position;
        }

        // -------- render --------
        window.clear(sf::Color(10, 10, 20));

        // trails
        if (trailsEnabled) {
            for (const auto& b : bodies) {
                if (b.trail.size() < 2) continue;
                std::vector<sf::Vertex> lineVertices;
                lineVertices.reserve(b.trail.size());
                for (const auto& p : b.trail) {
                    lineVertices.emplace_back(
                        worldToScreen(p, WIDTH, HEIGHT),
                        sf::Color(150, 150, 150)
                    );
                }
                window.draw(lineVertices.data(), lineVertices.size(), sf::LineStrip);
            }
        }

        // AddMoving arrow
        if (isDraggingAddMoving) {
            sf::Vector2f startScreen = worldToScreen(dragStart, WIDTH, HEIGHT);
            sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
            sf::Vertex line[] = {
                sf::Vertex(startScreen, sf::Color::White),
                sf::Vertex(mouseScreenF, sf::Color::White)
            };
            window.draw(line, 2, sf::Lines);
        }

        // bodies (zoom-aware radius & realistic colours)
        double zoomScale = 1.0e9 / metersPerPixel;
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            const Body& b = bodies[i];
            sf::Vector2f pos = worldToScreen(b.position, WIDTH, HEIGHT);

            float radius = (float)(b.radius * zoomScale);
            if (radius < 2.f)  radius = 2.f;
            if (radius > 60.f) radius = 60.f;

            sf::CircleShape c(radius);
            c.setOrigin(radius, radius);
            c.setPosition(pos);
            sf::Color col = bodyColor(b);
            c.setFillColor(col);

            if ((int)i == selectedIndex) {
                c.setOutlineThickness(2.f);
                c.setOutlineColor(sf::Color::Red);
            } else if (b.fixed) {
                c.setOutlineThickness(1.5f);
                c.setOutlineColor(sf::Color(230, 230, 80));
            } else if (b.ghost) {
                c.setOutlineThickness(1.5f);
                c.setOutlineColor(sf::Color(0, 150, 255));
            }
            window.draw(c);
        }

        // orbit overlay
        if (fontLoaded &&
            mode == InteractionMode::OrbitEdit &&
            orbitEdit.active &&
            orbitEdit.bodyIndex >= 0 &&
            orbitEdit.bodyIndex < (int)bodies.size() &&
            orbitEdit.parentIndex >= 0 &&
            orbitEdit.parentIndex < (int)bodies.size()) {
            drawOrbitOverlay(window,
                             bodies[orbitEdit.bodyIndex],
                             bodies[orbitEdit.parentIndex],
                             WIDTH, HEIGHT, font);
        }

        // HUD (top-left)
        if (fontLoaded) {
            sf::Text t;
            t.setFont(font);
            t.setCharacterSize(14);
            t.setFillColor(sf::Color::White);
            t.setPosition(10.f, 10.f);

            std::string modeStr;
            if (mode == InteractionMode::AddStill)       modeStr = "Add Still (1)";
            else if (mode == InteractionMode::AddMoving) modeStr = "Add Moving (2)";
            else                                         modeStr = "Add Orbiting (3)";

            std::string selStr = "None";
            if (selectedIndex >= 0 && selectedIndex < (int)bodies.size())
                selStr = bodies[selectedIndex].name;

            std::string hud = "Mode: " + modeStr +
                              "   Selected: " + selStr +
                              (lockVolume ? "   [Volume LOCKED]" : "   [Volume UNLOCKED]");
            t.setString(hud);
            window.draw(t);

            // bottom-left: scale & time
            sf::Text s1;
            s1.setFont(font);
            s1.setCharacterSize(14);
            s1.setFillColor(sf::Color::White);
            s1.setPosition(10.f, HEIGHT - 40.f);
            s1.setString("Scale: " + formatScale());
            window.draw(s1);

            sf::Text s2;
            s2.setFont(font);
            s2.setCharacterSize(14);
            s2.setFillColor(sf::Color::White);
            s2.setPosition(10.f, HEIGHT - 22.f);
            s2.setString("Time:  " + formatTimeScaleStr());
            window.draw(s2);
        }

        // Controls panel
        if (fontLoaded && controlsPanel.visible) {
            sf::RectangleShape panel;
            panel.setPosition(controlsPanel.pos);
            panel.setSize(controlsPanel.size);
            panel.setFillColor(sf::Color(20, 20, 35, 230));
            panel.setOutlineThickness(1.f);
            panel.setOutlineColor(sf::Color(120, 120, 120));
            window.draw(panel);

            float titleH = 24.f;
            sf::RectangleShape titleBar({controlsPanel.size.x, titleH});
            titleBar.setPosition(controlsPanel.pos);
            titleBar.setFillColor(sf::Color(40, 40, 60));
            window.draw(titleBar);

            sf::Text title;
            title.setFont(font);
            title.setCharacterSize(14);
            title.setFillColor(sf::Color::White);
            title.setPosition(controlsPanel.pos.x + 6.f, controlsPanel.pos.y + 4.f);
            title.setString("Controls  (? to toggle)");
            window.draw(title);

            float btnSize = 16.f;
            sf::RectangleShape minBtn({btnSize, btnSize});
            minBtn.setPosition(controlsPanel.pos.x + controlsPanel.size.x - (btnSize + 4.f),
                               controlsPanel.pos.y + 4.f);
            minBtn.setFillColor(sf::Color(80, 80, 80));
            window.draw(minBtn);

            sf::Text mText;
            mText.setFont(font);
            mText.setCharacterSize(12);
            mText.setFillColor(sf::Color::White);
            mText.setString(controlsPanel.minimized ? "+" : "-");
            mText.setPosition(minBtn.getPosition().x + 4.f, minBtn.getPosition().y - 3.f);
            window.draw(mText);

            if (!controlsPanel.minimized) {
                sf::Text line;
                line.setFont(font);
                line.setCharacterSize(12);
                line.setFillColor(sf::Color::White);
                float x = controlsPanel.pos.x + 8.f;
                float y = controlsPanel.pos.y + titleH + 6.f;

                auto l = [&](const std::string& s) {
                    line.setPosition(x, y);
                    line.setString(s);
                    window.draw(line);
                    y += 14.f;
                };

                l("LMB: add body (mode dependent)");
                l("RMB: select body (open editor)");
                l("1/2/3: Add Still / Moving / Orbiting");
                l("Enter (orbit edit): commit orbit");
                l("Esc: cancel rename/orbit/reset");
                l("Del: delete selected body");
                l("Up/Down: scale mass of selected");
                l("T: toggle trails   [ ]: trail length");
                l("V: toggle volume lock   X: toggle fixed");
                l("F: follow selected body   F2: rename");
                l("Ctrl+R: reset system (with confirm)");
                l("I/O: zoom in/out");
                l(", / . : slower / faster time");
            }
        }

        // Body editor panel
        if (fontLoaded && bodyPanel.visible &&
            selectedIndex >= 0 && selectedIndex < (int)bodies.size()) {
            Body& b = bodies[selectedIndex];

            sf::RectangleShape panel;
            panel.setPosition(bodyPanel.pos);
            panel.setSize(bodyPanel.size);
            panel.setFillColor(sf::Color(20, 20, 35, 230));
            panel.setOutlineThickness(1.f);
            panel.setOutlineColor(sf::Color::White);
            window.draw(panel);

            float titleH = 24.f;
            sf::RectangleShape titleBar({bodyPanel.size.x, titleH});
            titleBar.setPosition(bodyPanel.pos);
            titleBar.setFillColor(sf::Color(40, 40, 60));
            window.draw(titleBar);

            sf::Text title;
            title.setFont(font);
            title.setCharacterSize(14);
            title.setFillColor(sf::Color::White);
            title.setPosition(bodyPanel.pos.x + 6.f, bodyPanel.pos.y + 4.f);
            title.setString("Body Editor");
            window.draw(title);

            float btnSize = 16.f;
            sf::RectangleShape closeBtn({btnSize, btnSize});
            closeBtn.setPosition(bodyPanel.pos.x + bodyPanel.size.x - btnSize - 4.f,
                                 bodyPanel.pos.y + 4.f);
            closeBtn.setFillColor(sf::Color(120, 40, 40));
            window.draw(closeBtn);

            sf::RectangleShape minBtn({btnSize, btnSize});
            minBtn.setPosition(bodyPanel.pos.x + bodyPanel.size.x - 2.f * (btnSize + 4.f),
                               bodyPanel.pos.y + 4.f);
            minBtn.setFillColor(sf::Color(80, 80, 80));
            window.draw(minBtn);

            sf::Text xText;
            xText.setFont(font);
            xText.setCharacterSize(12);
            xText.setFillColor(sf::Color::White);
            xText.setString("x");
            xText.setPosition(closeBtn.getPosition().x + 4.f,
                              closeBtn.getPosition().y - 2.f);
            window.draw(xText);

            sf::Text mText;
            mText.setFont(font);
            mText.setCharacterSize(12);
            mText.setFillColor(sf::Color::White);
            mText.setString(bodyPanel.minimized ? "+" : "-");
            mText.setPosition(minBtn.getPosition().x + 4.f,
                              minBtn.getPosition().y - 3.f);
            window.draw(mText);

            if (!bodyPanel.minimized) {
                sf::Text line;
                line.setFont(font);
                line.setCharacterSize(13);
                line.setFillColor(sf::Color::White);
                float x = bodyPanel.pos.x + 8.f;
                float y = bodyPanel.pos.y + titleH + 6.f;

                auto drawLine = [&](const std::string& s) {
                    line.setPosition(x, y);
                    line.setString(s);
                    window.draw(line);
                    y += 20.f;
                };

                // name / type
                if (renaming)
                    drawLine("Name: " + renameBuffer + "_");
                else
                    drawLine("Name: " + b.name + "  (F2 to rename)");

                drawLine("Type: " + b.typeInfo.name);

                // simple type hint / cycling
                drawLine("Change type: [<] [>] (not clickable yet)");

                // mass (editable via +/-)
                drawLine("Mass: " + formatSci(b.mass) + " kg");

                // density & radius
                drawLine("Density: " + formatSci(b.density));
                drawLine("Base radius: " + formatSci(b.radius) + " px");

                // temperature
                if (isStar(b)) {
                    drawLine("Star T: " + formatSci(b.temperatureK) + " K");
                    drawLine("Lum: " + formatSci(b.luminosityW) + " W");
                } else {
                    drawLine("Surface T: " + formatSci(b.temperatureK) + " K");
                    drawLine("Albedo: " + formatSci(b.albedo));
                }

                // composition summary
                char buf[128];
                std::snprintf(buf, sizeof(buf),
                              "Comp R/M/I/G/W: %.2f/%.2f/%.2f/%.2f/%.2f",
                              b.comp.rock, b.comp.metal, b.comp.ice,
                              b.comp.gas, b.comp.water);
                drawLine(buf);

                // flags
                std::string flags = "Fixed: ";
                flags += (b.fixed ? "Yes" : "No");
                flags += "   Volume lock: ";
                flags += (lockVolume ? "On" : "Off");
                drawLine(flags);
            }
        }

        // Reset popup
        if (resetPopup.visible && fontLoaded) {
            float pw = 300.f, ph = 140.f;
            float px = (WIDTH - pw) * 0.5f;
            float py = (HEIGHT - ph) * 0.5f;

            sf::RectangleShape box({pw, ph});
            box.setPosition(px, py);
            box.setFillColor(sf::Color(20, 20, 35, 240));
            box.setOutlineThickness(1.f);
            box.setOutlineColor(sf::Color::White);
            window.draw(box);

            sf::Text msg;
            msg.setFont(font);
            msg.setCharacterSize(16);
            msg.setFillColor(sf::Color::White);
            msg.setString("Reset system to only the star?");
            msg.setPosition(px + 20.f, py + 20.f);
            window.draw(msg);

            float btnW = 80.f, btnH = 28.f, gap = 20.f;
            float yesX = px + pw * 0.5f - gap - btnW;
            float noX  = px + pw * 0.5f + gap;
            float btnY = py + ph - btnH - 20.f;

            sf::RectangleShape yesBtn({btnW, btnH});
            yesBtn.setPosition(yesX, btnY);
            yesBtn.setFillColor(sf::Color(60, 120, 60));
            window.draw(yesBtn);

            sf::RectangleShape noBtn({btnW, btnH});
            noBtn.setPosition(noX, btnY);
            noBtn.setFillColor(sf::Color(120, 60, 60));
            window.draw(noBtn);

            sf::Text yesTxt;
            yesTxt.setFont(font);
            yesTxt.setCharacterSize(14);
            yesTxt.setFillColor(sf::Color::White);
            yesTxt.setString("Yes");
            yesTxt.setPosition(yesX + 24.f, btnY + 4.f);
            window.draw(yesTxt);

            sf::Text noTxt;
            noTxt.setFont(font);
            noTxt.setCharacterSize(14);
            noTxt.setFillColor(sf::Color::White);
            noTxt.setString("No");
            noTxt.setPosition(noX + 28.f, btnY + 4.f);
            window.draw(noTxt);

            resetPopup.yesRect = sf::FloatRect(yesX, btnY, btnW, btnH);
            resetPopup.noRect  = sf::FloatRect(noX,  btnY, btnW, btnH);
        }

        window.display();
    }

    return 0;
}
