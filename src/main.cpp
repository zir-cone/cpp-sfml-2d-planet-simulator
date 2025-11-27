#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <cstdio>
#include <algorithm>

// --------------------
// Basic 2D vector
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
// Im Gonna Kill Myself <3 <3 <3 xoxoxo

// --------------------
// Types & composition
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

struct CompositionBulk {
    double rock  = 0.0;
    double metal = 0.0;
    double ice   = 0.0;
    double gas   = 0.0;
    double water = 0.0;
};

struct SubComponent {
    std::string name;
    double fraction; // mass fraction (0..1)
};

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
// Body
// --------------------
struct Body {
    double mass;         // kg
    Vec2   position;     // m
    Vec2   velocity;     // m/s
    Vec2   acceleration; // m/s^2

    double density;
    double radius;       // base visual radius at reference zoom (px)
    bool   fixed  = false;
    bool   ghost  = false;

    std::string  name;
    BodyTypeInfo typeInfo;
    BodyCategory category = BodyCategory::Rocky;

    CompositionBulk bulk;
    std::vector<SubComponent> sub;

    // thermal / radiative
    double temperatureK = 288.0;   // surface temp (planets) or effective temp (stars)
    double luminosityW  = 0.0;     // stars only
    double albedo       = 0.3;     // planets

    sf::Color color;
    std::vector<Vec2> trail;
};

// --------------------
// Global constants
// --------------------
const double G_PHYS       = 6.67430e-11;              // m^3 kg^-1 s^-2
const double SIGMA_SB     = 5.670374419e-8;           // W m^-2 K^-4
const double YEAR_SECONDS = 365.25 * 24.0 * 3600.0;   // Yuh, what is this diddy blud doing on the calculator, is blud Einsten?

const double M_SUN   = 1.98847e30;
const double M_EARTH = 5.97219e24;

// scale / time (mutable)
double metersPerPixel = 1.0e9;                        // 1 px = 1e9 m
double timeScale      = 86400.0 * 50.0;               // 50 days per real second

const double SOFTENING_METERS = 1.0e9;

bool trailsEnabled = true;
int  maxTrailPoints = 6700;                           // SIX SEVENENNNE NUINomjiknF OW
bool lockVolume     = false;                          // ive lost my shit

// helpers radius/density
double computeVolumeFromRadius(double r) {
    const double pi = 3.141592653589793;            // What is this diddy blud doing on the calculator? Is blud Einsten? What is this diddyblud doing on the calculator? Does he think he is, Epstein?
    return 4.0 / 3.0 * pi * r * r * r;
}
double computeRadiusFromMassDensity(double mass, double density) {
    if (density <= 0.0) density = 1.0;
    double volume = mass / density;
    const double pi = 3.141592653589793;           // pi lwk tuff bro lk tht 3.14159 vro ra 🌹💖
    return std::cbrt((3.0 * volume) / (4.0 * pi));
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
    return sf::Vector2f(static_cast<float>(cx + dx),
                        static_cast<float>(cy + dy));
}

// --------------------
// Gravity & thermal
// --------------------
bool isStar(const Body& b) {
    return b.category == BodyCategory::Star;
}

void computeGravity(std::vector<Body>& bodies) {
    std::size_t n = bodies.size();
    for (auto& b : bodies) b.acceleration = {0.0, 0.0};

    for (std::size_t i = 0; i < n; ++i) {
        if (bodies[i].ghost) continue;
        for (std::size_t j = i + 1; j < n; ++j) {
            if (bodies[j].ghost) continue;
            Vec2 r = bodies[j].position - bodies[i].position;
            double dist2 = r.x * r.x + r.y * r.y;
            double soft2 = dist2 + SOFTENING_METERS * SOFTENING_METERS;
            double soft  = std::sqrt(soft2);
            if (soft < 1.0) soft = 1.0;
            double invDist3 = 1.0 / (soft2 * soft);

            Vec2 acc_i = r * (G_PHYS * bodies[j].mass * invDist3);
            Vec2 acc_j = r * (-G_PHYS * bodies[i].mass * invDist3);
            bodies[i].acceleration += acc_i;
            bodies[j].acceleration += acc_j;
        }
    }
}

void updateTemperatures(std::vector<Body>& bodies) {
    std::vector<int> starIdx;
    for (int i = 0; i < (int)bodies.size(); ++i)
        if (isStar(bodies[i]) && bodies[i].luminosityW > 0.0)
            starIdx.push_back(i);
    if (starIdx.empty()) return;

    for (int i = 0; i < (int)bodies.size(); ++i) {
        Body& body = bodies[i];
        if (isStar(body)) continue;

        double bestDist2 = 1e99;
        const Body* bestStar = nullptr;
        for (int si : starIdx) {
            const Body& s = bodies[si];
            Vec2 r = body.position - s.position;
            double d2 = r.x*r.x + r.y*r.y;
            if (d2 < bestDist2) { bestDist2 = d2; bestStar = &s; }
        }
        if (!bestStar) continue;

        double r = std::sqrt(bestDist2);
        if (r <= 0.0) r = 1.0;
        double F = bestStar->luminosityW / (4.0 * 3.141592653589793 * r * r);
        double T = std::pow((F * (1.0 - body.albedo)) / (4.0 * SIGMA_SB), 0.25);
        if (std::isfinite(T)) body.temperatureK = T;
    }
}

int findNearestStar(const std::vector<Body>& bodies, int bodyIndex) {
    int best = -1;
    double bestDist2 = 1e99;
    for (int i = 0; i < (int)bodies.size(); ++i) {
        if (i == bodyIndex) continue;
        if (!isStar(bodies[i])) continue;
        Vec2 d{ bodies[i].position.x - bodies[bodyIndex].position.x,
                bodies[i].position.y - bodies[bodyIndex].position.y };
        double d2 = d.x*d.x + d.y*d.y;
        if (d2 < bestDist2) {
            bestDist2 = d2;
            best = i;
        }
    }
    return best;
}

// --------------------
// Star & planet category helpers
// --------------------
void updateStarFromMass(Body& b) {
    if (!isStar(b)) return;
    double mMsun = b.mass / M_SUN;
    if (mMsun < 0.1) mMsun = 0.1;
    if (mMsun > 50.0) mMsun = 50.0;

    const double Lsun = 3.828e26;
    const double Tsun = 5800.0;

    b.luminosityW  = Lsun * std::pow(mMsun, 3.5);
    b.temperatureK = Tsun * std::pow(mMsun, 0.5);

    double Rsun_px = 18.0;
    b.radius = Rsun_px * std::pow(mMsun, 0.8);
}

sf::Color colorFromStarTemp(double T) {
    T = std::clamp(T, 2000.0, 40000.0);
    double t = T / 100.0;
    double r, g, b;
    if (t <= 66.0) {
        r = 255.0;
        g = 99.4708025861 * std::log(t) - 161.1195681661;
        if (t <= 19.0) b = 0.0;
        else          b = 138.5177312231 * std::log(t - 10.0) - 305.0447927307;
    } else {
        r = 329.698727446 * std::pow(t - 60.0, -0.1332047592);
        g = 288.1221695283 * std::pow(t - 60.0, -0.0755148492);
        b = 255.0;
    }
    auto cl = [](double x) {
        if (x < 0.0) x = 0.0;
        if (x > 255.0) x = 255.0;
        return (sf::Uint8)x;
    };
    return sf::Color(cl(r), cl(g), cl(b));
}

sf::Color colorFromPlanet(const Body& b) {
    sf::Color base;

    double water = b.bulk.water;
    double ice   = b.bulk.ice;
    double gas   = b.bulk.gas;

    if (b.category == BodyCategory::GasGiant) {
        // tan / reddish gas-giant palette
        base = sf::Color(205, 175, 130);
    } else if (b.category == BodyCategory::IceGiant) {
        // pale blue
        base = sf::Color(115, 175, 230);
    } else {
        // Rock-likes: water-heavy => oceanic blue, else rocky / metallic browns
        if (water > 0.25)
            base = sf::Color(60, 110, 170);      // marine / oceanic
        else if (ice > 0.3)
            base = sf::Color(180, 220, 250);     // icy silicate
        else
            base = sf::Color(150, 120, 90);      // generic rocky
    }

    // simple temperature tinting
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
    return isStar(b) ? colorFromStarTemp(b.temperatureK)
                     : colorFromPlanet(b);
}

void updatePlanetCategoryFromComp(Body& b) {
    if (isStar(b)) return;

    double rock   = b.bulk.rock;
    double metal  = b.bulk.metal;
    double ice    = b.bulk.ice;
    double gas    = b.bulk.gas;
    double water  = b.bulk.water;

    double mEarth = b.mass / M_EARTH;
    double mJup   = b.mass / 1.898e27; // Jupiter mass

    // size label
    std::string sizeLabel;
    if (mEarth < 0.5)       sizeLabel = "subterra";
    else if (mEarth < 2.0)  sizeLabel = "terra";
    else if (mEarth < 10.0) sizeLabel = "superterra";
    else                    sizeLabel = "superearth";

    bool nontectonic = (mEarth < 0.4);   // crude tectonics heuristic
    bool marine      = (water > 0.25 && (rock + metal) > 0.4);
    bool metallic    = (metal > rock);

    // --- Gas / ice giant side ---
    if (gas > 0.6 && mJup >= 0.1) {
        b.category = BodyCategory::GasGiant;

        if (mJup < 0.5)       b.typeInfo.name = "Gaseous subjupiter";
        else if (mJup < 2.0)  b.typeInfo.name = "Jovian";
        else                  b.typeInfo.name = "Superjupiter";

        b.albedo = 0.4;
    }
    else if ((ice + water + gas) > 0.6 && mEarth > 5.0) {
        b.category = BodyCategory::IceGiant;

        if (mEarth < 10.0)       b.typeInfo.name = "Ice giant subneptune";
        else if (mEarth < 25.0)  b.typeInfo.name = "Ice giant";
        else                     b.typeInfo.name = "Superneptune";

        b.albedo = 0.5;
    }
    // --- Terrestrial side ---
    else {
        b.category = BodyCategory::Rocky;

        std::string prefix;
        if (marine)          prefix = "Marine continental ";
        else if (metallic)   prefix = "Metallic ";
        else                 prefix = "Rocky ";

        std::string tect;
        if (nontectonic) tect = "nontectonic ";

        b.typeInfo.name = prefix + tect + sizeLabel;
        b.albedo = marine ? 0.3 : 0.25;
    }

    b.typeInfo.description = "Auto-classified from mass and bulk composition.";
    b.color = bodyColor(b);
}

// Map detailed subcomponents into bulk categories
void recomputeBulkComposition(Body& b) {
    b.bulk = {};
    double sum = 0.0;
    for (const auto& sc : b.sub) sum += sc.fraction;
    if (sum <= 0.0) return;

    for (const auto& sc : b.sub) {
        double f = sc.fraction / sum;
        std::string n = sc.name;
        if (n == "Silicate rock")         b.bulk.rock  += f;
        else if (n == "Iron")             b.bulk.metal += f;
        else if (n == "Water")            b.bulk.water += f;
        else if (n == "Ice")              b.bulk.ice   += f;
        else if (n == "Hydrogen" ||
                 n == "Helium"   ||
                 n == "Nitrogen" ||
                 n == "Oxygen"   ||
                 n == "CO2"      ||
                 n == "Methane"  ||
                 n == "Argon"    ||
                 n == "Metals")           b.bulk.gas   += f;
        else                               b.bulk.rock  += f * 0.5;
    }

    if (!isStar(b))
        updatePlanetCategoryFromComp(b);
}

// --------------------
// Integrator
// --------------------
// const double SECONDS_PER_DAY = 86400.0; im a retard for multiplying this with timescale

void step(std::vector<Body>& bodies, double dtRealSeconds)
{
    // total simulated time to advance this frame
    double dtSim = dtRealSeconds * timeScale;

    // maximum allowed step
    const double MAX_SUBSTEP = 3600.0; // 1 simulated hour per substep

    double absDt = std::abs(dtSim);
    int nSub = (absDt <= MAX_SUBSTEP) ? 1
                                      : (int)std::ceil(absDt / MAX_SUBSTEP); // i dont know why i do shit like this at 4 am, like why the FUCK is it split like this????
                                                                             // im gonna keep it though cause its quirky and the cup runeth over on quirk, wow!

    double h = dtSim / nSub;  // size of each substep

    bool recordTrails = trailsEnabled;
    if (!recordTrails) {
        for (auto& b : bodies)
            b.trail.clear();
    }

    for (int s = 0; s < nSub; ++s) {
        // v_{n+1} = v_n + a_n * h
        for (auto& b : bodies) {
            if (!b.fixed && !b.ghost) {
                b.velocity += b.acceleration * h;
            }
        }

        // x_{n+1} = x_n + v_{n+1} * h  (semi-implicit Euler)
        for (auto& b : bodies) {
            if (!b.fixed && !b.ghost) {
                b.position += b.velocity * h;
            }
        }

        // update accelerations for next substep
        computeGravity(bodies);

        if (recordTrails) {
            for (auto& b : bodies) {
                b.trail.push_back(b.position);
                if ((int)b.trail.size() > maxTrailPoints) {
                    b.trail.erase(b.trail.begin());
                }
            }
        }
    }
}

// --------------------
// Parent selection
// --------------------
int findStrongestInfluenceParent(const std::vector<Body>& bodies,
                                 const Vec2& posMeters)
{
    int bestIndex = -1;
    double bestAccel = -1.0;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        const Body& b = bodies[i];
        if (b.ghost) continue;
        Vec2 r = b.position - posMeters;
        double dist2 = r.x * r.x + r.y * r.y;
        if (dist2 <= 0.0) continue;
        double accel = G_PHYS * b.mass / dist2;
        if (accel > bestAccel) {
            bestAccel = accel;
            bestIndex = (int)i;
        }
    }
    return bestIndex;
}

// nvm ts back to a shitfuck 😭
int chooseOrbitParent(const std::vector<Body>& bodies, const Vec2& posMeters) {
    int bestPlanet = -1;
    double bestScore = 1e99;

    // first pass
    for (int i = 0; i < (int)bodies.size(); ++i) {
        const Body& p = bodies[i];
        if (isStar(p)) continue;

        // nearest star for this planet
        int starIdx = findNearestStar(bodies, i);
        if (starIdx < 0) continue;
        const Body& star = bodies[starIdx];

        // planet-star distance
        Vec2 ps = p.position - star.position;
        double dStar = length(ps);
        if (dStar <= 0.0) continue;

        // Hill radius
        double rHill = dStar * std::cbrt(p.mass / (3.0 * star.mass));
        if (rHill <= 0.0) continue;

        // distance from click to planet
        Vec2 dp = posMeters - p.position;
        double dMoon = length(dp);

        if (dMoon < rHill) {
            // "i can feel it deep inside deep deep down inside" - Fade by Kanye West (The Life of Pablo, 2016)
            double score = dMoon / rHill;
            if (score < bestScore) {
                bestScore = score;
                bestPlanet = i;
            }
        }
    }

    if (bestPlanet != -1) {
        return bestPlanet;  // click is inside some planet’s Hill sphere → moon
    }

    // Otherwise: pick closest star as parent
    int bestStar = -1;
    double bestD2 = 1e99;
    for (int i = 0; i < (int)bodies.size(); ++i) {
        if (!isStar(bodies[i])) continue;
        Vec2 ds = posMeters - bodies[i].position;
        double d2 = ds.x*ds.x + ds.y*ds.y;
        if (d2 < bestD2) {
            bestD2 = d2;
            bestStar = i;
        }
    }

    if (bestStar != -1)
        return bestStar;

    // fallback
    return findStrongestInfluenceParent(bodies, posMeters);
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

int findNearestBody(const std::vector<Body>& bodies,
                    const Vec2& pointMeters,
                    double maxDistMeters) {
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

// --------------------
// Panels
// --------------------
struct Panel {
    bool visible = false;
    bool minimized = false;
    sf::Vector2f pos{20.0f, 80.0f};
    sf::Vector2f size{360.0f, 280.0f};
    bool dragging = false;
    sf::Vector2f dragOffset{0.0f, 0.0f};
    float scroll = 0.0f;   // vertical scroll offset for content
};

bool pointInRect(const sf::Vector2f& p, const sf::FloatRect& r) {
    return r.contains(p);
}

// --------------------
// Orbit edit
// --------------------
enum class OrbitHandle { None, Body, Pe, Ap };

struct OrbitEditState {
    bool   active = false;
    int    bodyIndex   = -1;
    int    parentIndex = -1;
    double a  = 0.0;
    double e  = 0.0;
    double nu = 0.0;
    Vec2   u{1.0, 0.0};
    Vec2   v{0.0, 1.0};
    OrbitHandle dragHandle = OrbitHandle::None;
    bool   dragging = false;
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
    body.position = { parent.position.x + rel.x,
                      parent.position.y + rel.y };

    double mu = G_PHYS * (body.mass + parent.mass);
    double vmag = std::sqrt(std::max(0.0, mu * (2.0 / r - 1.0 / a)));
    Vec2 rhat = rel; normalize(rhat);
    Vec2 that{ -rhat.y, rhat.x };

    body.velocity = { parent.velocity.x + that.x * vmag,
                      parent.velocity.y + that.y * vmag };
}

void initOrbitEdit(int bodyIdx, int parentIdx, std::vector<Body>& bodies) {
    orbitEdit.active = true;
    orbitEdit.bodyIndex = bodyIdx;
    orbitEdit.parentIndex = parentIdx;
    orbitEdit.dragging = false;
    orbitEdit.dragHandle = OrbitHandle::None;

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

OrbitHandle pickOrbitHandle(const Body& body,
                            const Body& parent,
                            const sf::Vector2f& mouseScreen,
                            double width,
                            double height)
{
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
        if (orbitEdit.e < 0.0)  orbitEdit.e = 0.0;
        if (orbitEdit.e > 0.99) orbitEdit.e = 0.99;

        orbitEdit.nu = angleMouse;
    }

    updateBodyFromOrbit(body, parent);
}
// Assassinate Yakub
// --------------------
// Draw helpers
// --------------------
void drawDashedLine(sf::RenderWindow& window,
                    const sf::Vector2f& A,
                    const sf::Vector2f& B,
                    sf::Color col = sf::Color::White) {
    const float dashLen = 8.f;
    const float gapLen  = 4.f;
    sf::Vector2f d = B - A;
    float len = std::sqrt(d.x*d.x + d.y*d.y);
    if (len <= 0.f) return;
    sf::Vector2f dir = d / len;
    float pos = 0.f;
    while (pos < len) {
        float seg = std::min(dashLen, len - pos);
        sf::Vector2f s = A + dir * pos;
        sf::Vector2f e = A + dir * (pos + seg);
        sf::Vertex line[] = { sf::Vertex(s, col), sf::Vertex(e, col) };
        window.draw(line, 2, sf::Lines);
        pos += dashLen + gapLen;
    }
}

// orbit overlay
void drawOrbitOverlay(sf::RenderWindow& window,
                      const Body& body,
                      const Body& parent,
                      double WIDTH,
                      double HEIGHT,
                      const sf::Font& font) {
    if (!orbitEdit.active) return;

    OrbitGeom geom;
    if (!computeOrbitGeom(parent, geom)) return;

    const int SEG = 256;
    std::vector<sf::Vertex> line;
    line.reserve(SEG + 1);
    for (int i = 0; i <= SEG; ++i) {
        double th = (double)i / SEG * 2.0 * 3.141592653589793;
        double cosT = std::cos(th);
        double sinT = std::sin(th);
        double denom = 1.0 + geom.e * cosT;
        if (denom <= 0.0) denom = 1e-6;
        double r = geom.a * (1.0 - geom.e*geom.e) / denom;
        Vec2 rel = geom.u * (r * cosT) + geom.v * (r * sinT);
        Vec2 pWorld{ geom.parentPos.x + rel.x,
                     geom.parentPos.y + rel.y };
        line.emplace_back(worldToScreen(pWorld, WIDTH, HEIGHT),
                          sf::Color(0, 150, 255));
    }
    window.draw(line.data(), line.size(), sf::LineStrip);

    sf::Vector2f parentScreen = worldToScreen(geom.parentPos, WIDTH, HEIGHT);
    sf::CircleShape parentPt(4.f);
    parentPt.setOrigin(4.f, 4.f);
    parentPt.setPosition(parentScreen);
    parentPt.setFillColor(sf::Color::Red);
    window.draw(parentPt);

    sf::Vector2f peScreen = worldToScreen(geom.peWorld, WIDTH, HEIGHT);
    sf::Vector2f apScreen = worldToScreen(geom.apWorld, WIDTH, HEIGHT);

    sf::CircleShape pePt(4.f), apPt(4.f);
    pePt.setOrigin(4.f, 4.f);
    apPt.setOrigin(4.f, 4.f);
    pePt.setPosition(peScreen);
    apPt.setPosition(apScreen);
    pePt.setFillColor(sf::Color(255, 100, 200));
    apPt.setFillColor(sf::Color(255, 100, 200));
    window.draw(pePt);
    window.draw(apPt);

    drawDashedLine(window, parentScreen, peScreen);
    drawDashedLine(window, parentScreen, apScreen);

    if (!font.getInfo().family.empty()) {
        sf::Text text;
        text.setFont(font);
        text.setCharacterSize(12);

        text.setFillColor(sf::Color(255, 100, 200));
        text.setString("Pe");
        text.setPosition(peScreen.x + 6.f, peScreen.y - 4.f);
        window.draw(text);
        text.setString("Ap");
        text.setPosition(apScreen.x + 6.f, apScreen.y - 4.f);
        window.draw(text);

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
        std::snprintf(buf, sizeof(buf), "nu = %.2f°", nuDeg);
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
void resetSystem(std::vector<Body>& bodies, Body& starProto) {
    bodies.clear();
    Body star = starProto;
    star.position = {0.0, 0.0};
    star.velocity = {0.0, 0.0};
    star.fixed    = true;
    star.ghost    = false;
    updateStarFromMass(star);
    star.color = bodyColor(star);
    bodies.push_back(star);
    computeGravity(bodies);
    updateTemperatures(bodies);
    for (auto& b : bodies) b.trail.clear();
    orbitEdit = {};
    viewBodyIndex = -1;
    viewCenterWorld = star.position;
}

// --------------------
// Body panel controls / rects
// --------------------
struct BodyPanelControls {
    bool valid = false;
    sf::FloatRect compRock;
    sf::FloatRect compMetal;
    sf::FloatRect compIce;
    sf::FloatRect compGas;
    sf::FloatRect compWater;
};
BodyPanelControls bodyPanelControls;
std::vector<sf::FloatRect> subCompRects; // detailed composition rows

void setPrimaryComposition(Body& b, const std::string& which) {
    if (isStar(b)) {
        if (which == "Gas") {
            b.bulk = {0,0,0,1.0,0};
        } else if (which == "Metal") {
            b.bulk = {0,1.0,0,0,0};
        }
        updatePlanetCategoryFromComp(b);
        return;
    }

    b.bulk = {};
    if (which == "Rock")      b.bulk.rock  = 1.0;
    else if (which == "Metal") b.bulk.metal = 1.0;
    else if (which == "Ice")   b.bulk.ice   = 1.0;
    else if (which == "Gas")   b.bulk.gas   = 1.0;
    else if (which == "Water") b.bulk.water = 1.0;

    updatePlanetCategoryFromComp(b);
}

// -------- Lighting (shader) --------
sf::Texture gPlanetMask;
sf::Shader  gPlanetShader;
bool gHasLightingShader = false;

const char* PLANET_FRAGMENT_SHADER = R"(
uniform sampler2D texture;
uniform vec3 baseColor;
uniform vec2 lightDir;  // direction *to* the star
uniform float ambient;

void main()
{
    vec2 uv = gl_TexCoord[0].xy;
    vec4 mask = texture2D(texture, uv);
    if (mask.a < 0.05)
        discard;

    // map UV -> [-1,1] disk
    vec2 p = uv * 2.0 - 1.0;
    float r2 = dot(p, p);
    if (r2 > 1.0)
        discard;

    // reconstruct a 3D unit-sphere normal
    // note: screen Y goes down, so flip sign to get math-style up
    float z = sqrt(max(0.0, 1.0 - r2));
    vec3 n = normalize(vec3(p.x, -p.y, z));

    // light direction in the same space (z = 0)
    vec3 L = normalize(vec3(lightDir.x, -lightDir.y, 0.0));

    float ndotl = max(dot(n, L), 0.0);
    float lighting = ambient + (1.0 - ambient) * ndotl;

    vec3 col = baseColor * lighting;
    gl_FragColor = vec4(col, mask.a);
}
)";

// --------------------
// fml
// --------------------
int main() {
    const unsigned int WIDTH = 1280;
    const unsigned int HEIGHT = 720;

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "PlanetSim.exe");
    window.setFramerateLimit(120);

        // --- Lighting init ---
    if (sf::Shader::isAvailable()) {
        // circular alpha mask
        sf::Image img;
        const int TEX_SIZE = 256;
        img.create(TEX_SIZE, TEX_SIZE, sf::Color::Transparent);
        for (int y = 0; y < TEX_SIZE; ++y) {
            for (int x = 0; x < TEX_SIZE; ++x) {
                float fx = (x + 0.5f) / TEX_SIZE * 2.f - 1.f;
                float fy = (y + 0.5f) / TEX_SIZE * 2.f - 1.f;
                float r2 = fx*fx + fy*fy;
                if (r2 <= 1.f) {
                    img.setPixel(x, y, sf::Color(255, 255, 255, 255));
                }
            }
        }
        gPlanetMask.loadFromImage(img);
        gPlanetMask.setSmooth(true);

        gHasLightingShader =
            gPlanetShader.loadFromMemory(PLANET_FRAGMENT_SHADER, sf::Shader::Fragment);
    } else {
        gHasLightingShader = false;
    }

    std::vector<Body> bodies;

    // prototype star (for reset)
    Body starProto;
    starProto.mass        = M_SUN;
    starProto.position    = {0.0, 0.0};
    starProto.velocity    = {0.0, 0.0};
    starProto.fixed       = true;
    starProto.density     = 1.0;
    starProto.radius      = 18.0;
    starProto.name        = "Sol";
    starProto.category    = BodyCategory::Star;
    starProto.typeInfo    = starType;
    starProto.temperatureK = 5800.0;
    starProto.luminosityW  = 3.828e26;
    starProto.albedo       = 0.0;

    starProto.sub = {
        {"Hydrogen", 0.70},
        {"Helium",   0.28},
        {"Metals",   0.02}
    };
    recomputeBulkComposition(starProto);
    updateStarFromMass(starProto);
    starProto.color = bodyColor(starProto);

    bodies.push_back(starProto);
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
    bodyPanel.size = {360.f, 280.f};

    Panel controlsPanel;
    controlsPanel.pos  = {WIDTH - 380.f, 10.f};
    controlsPanel.size = {360.f, 260.f};
    controlsPanel.visible = true;

    bool renaming = false;
    std::string renameBuffer;

    viewCenterWorld = starProto.position;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {

            if (event.type == sf::Event::Closed) {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::R && event.key.control) {
                    resetPopup.visible = true;
                }
                else if (event.key.code == sf::Keyboard::Slash && event.key.shift) {
                    controlsPanel.visible = !controlsPanel.visible;
                }
                else if (event.key.code == sf::Keyboard::Escape) {
                    if (resetPopup.visible) {
                        resetPopup.visible = false;
                    } else if (renaming) {
                        renaming = false;
                    } else if (mode == InteractionMode::OrbitEdit && orbitEdit.active &&
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
                            Body& b = bodies[selectedIndex];
                            b.mass *= 1.2;
                            updateBodyAfterMassChange(b);
                            if (isStar(b)) {
                                updateStarFromMass(b);
                                b.color = bodyColor(b);
                            }
                        }
                    }
                    else if (event.key.code == sf::Keyboard::Down) {
                        if (selectedIndex >= 0) {
                            Body& b = bodies[selectedIndex];
                            b.mass *= 0.8;
                            updateBodyAfterMassChange(b);
                            if (isStar(b)) {
                                updateStarFromMass(b);
                                b.color = bodyColor(b);
                            }
                        }
                    }
                    else if (event.key.code == sf::Keyboard::T) {
                        trailsEnabled = !trailsEnabled;
                    }
                    else if (event.key.code == sf::Keyboard::LBracket) {
                        if (maxTrailPoints > 100) { maxTrailPoints -= 100;} 
                    }
                    else if (event.key.code == sf::Keyboard::RBracket) {
                        maxTrailPoints += 100;
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
            else if (event.type == sf::Event::MouseWheelScrolled) {
                sf::Vector2f mouse = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                sf::FloatRect panelRect(bodyPanel.pos, bodyPanel.size);

                if (panelRect.contains(mouse) && bodyPanel.visible && !bodyPanel.minimized) {
                    // delta > 0 => scroll up
                    bodyPanel.scroll += event.mouseWheelScroll.delta * 20.f;

                    // clamp scroll so we don't wander off into nowhere
                    const float minScroll = -260.f; // tune as needed
                    const float maxScroll = 0.f;
                    if (bodyPanel.scroll > maxScroll) bodyPanel.scroll = maxScroll;
                    if (bodyPanel.scroll < minScroll) bodyPanel.scroll = minScroll;
                }
            }
            else if (event.type == sf::Event::MouseButtonPressed) {
                sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));

                // reset popup
                if (resetPopup.visible && event.mouseButton.button == sf::Mouse::Left) {
                    if (resetPopup.yesRect.contains(mouseScreenF)) {
                        resetPopup.visible = false;
                        resetSystem(bodies, starProto);
                    } else if (resetPopup.noRect.contains(mouseScreenF)) {
                        resetPopup.visible = false;
                    }
                    continue;
                }

                // ----- Body panel handling -----
                if (bodyPanel.visible) {
                    sf::FloatRect panelRect(bodyPanel.pos, bodyPanel.size);
                    sf::FloatRect titleRect(bodyPanel.pos.x, bodyPanel.pos.y,
                                            bodyPanel.size.x, 24.f);

                    if (event.mouseButton.button == sf::Mouse::Left ||
                        event.mouseButton.button == sf::Mouse::Right) {
                        if (titleRect.contains(mouseScreenF)) {
                            float btn = 16.f;
                            sf::FloatRect closeRect(
                                bodyPanel.pos.x + bodyPanel.size.x - btn - 4.f,
                                bodyPanel.pos.y + 4.f, btn, btn);
                            sf::FloatRect minRect(
                                bodyPanel.pos.x + bodyPanel.size.x - 2.f*(btn+4.f),
                                bodyPanel.pos.y + 4.f, btn, btn);
                            if (closeRect.contains(mouseScreenF) &&
                                event.mouseButton.button == sf::Mouse::Left) {
                                bodyPanel.visible = false;
                                selectedIndex = -1;
                                continue;
                            } else if (minRect.contains(mouseScreenF) &&
                                       event.mouseButton.button == sf::Mouse::Left) {
                                bodyPanel.minimized = !bodyPanel.minimized;
                                continue;
                            } else if (event.mouseButton.button == sf::Mouse::Left) {
                                bodyPanel.dragging = true;
                                bodyPanel.dragOffset = mouseScreenF - bodyPanel.pos;
                                continue;
                            }
                            // if somebody other than me ever happens to read this code
                            // and doesnt know c++ please never learn this language
                            // if you do know c++ please unlearn this language
                            // this shit the type of pain used in hell 🥀🥀🥀
                        } else if (panelRect.contains(mouseScreenF)) {
                            // content click: primary comp boxes
                            if (selectedIndex >= 0 &&
                                selectedIndex < (int)bodies.size() &&
                                bodyPanelControls.valid &&
                                event.mouseButton.button == sf::Mouse::Left) {

                                Body& b = bodies[selectedIndex];
                                if (bodyPanelControls.compRock.contains(mouseScreenF)) {
                                    setPrimaryComposition(b, "Rock");
                                } else if (bodyPanelControls.compMetal.contains(mouseScreenF)) {
                                    setPrimaryComposition(b, "Metal");
                                } else if (bodyPanelControls.compIce.contains(mouseScreenF)) {
                                    setPrimaryComposition(b, "Ice");
                                } else if (bodyPanelControls.compGas.contains(mouseScreenF)) {
                                    setPrimaryComposition(b, "Gas");
                                } else if (bodyPanelControls.compWater.contains(mouseScreenF)) {
                                    setPrimaryComposition(b, "Water");
                                }
                            }

                            // detailed subcomponents: left click +0.05, right click -0.05
                            if (selectedIndex >= 0 &&
                                selectedIndex < (int)bodies.size() &&
                                !subCompRects.empty()) {

                                Body& b = bodies[selectedIndex];
                                for (std::size_t i = 0; i < subCompRects.size(); ++i) {
                                    if (subCompRects[i].contains(mouseScreenF)) {
                                        if (i >= b.sub.size()) break;
                                        double delta = (event.mouseButton.button == sf::Mouse::Left)
                                                       ? +0.05 : -0.05;
                                        b.sub[i].fraction = std::clamp(b.sub[i].fraction + delta, 0.0, 1.0);
                                        // normalize
                                        double sum = 0.0;
                                        for (auto& sc : b.sub) sum += sc.fraction;
                                        if (sum > 0.0) {
                                            for (auto& sc : b.sub) sc.fraction /= sum;
                                        }
                                        recomputeBulkComposition(b);
                                        b.color = bodyColor(b);
                                        break;
                                    }
                                }
                            }
                            continue;
                        }
                    }
                }
                // ----- Controls panel -----
                if (controlsPanel.visible) {
                    sf::FloatRect titleRect(controlsPanel.pos.x, controlsPanel.pos.y,
                                            controlsPanel.size.x, 24.f);
                    sf::FloatRect panelRect(controlsPanel.pos, controlsPanel.size);

                    if (event.mouseButton.button == sf::Mouse::Left) {
                        if (titleRect.contains(mouseScreenF)) {
                            float btn = 16.f;
                            sf::FloatRect minRect(
                                controlsPanel.pos.x + controlsPanel.size.x - (btn + 4.f),
                                controlsPanel.pos.y + 4.f, btn, btn);
                            if (minRect.contains(mouseScreenF)) {
                                controlsPanel.minimized = !controlsPanel.minimized;
                                continue;
                            } else {
                                controlsPanel.dragging = true;
                                controlsPanel.dragOffset = mouseScreenF - controlsPanel.pos;
                                continue;
                            }
                        } else if (panelRect.contains(mouseScreenF)) {
                            continue; // absorb
                        }
                    }
                }

                // ----- Orbit handle pick -----
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

                // ----- Main mouse logic -----
                double cx = WIDTH * 0.5;
                double cy = HEIGHT * 0.5;
                Vec2 worldPos{
                    (mouseScreenF.x - cx) * metersPerPixel + viewCenterWorld.x,
                    (mouseScreenF.y - cy) * metersPerPixel + viewCenterWorld.y
                };

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
                        b.albedo   = 0.25;
                        b.name     = "Rocky body";
                        b.sub = {
                            {"Silicate rock", 0.67}, // 67 🥀🥀🥀🥀🥀🥀🥀🥀🥀🥀
                            {"Iron",          0.32},
                            {"Water",         0.01},
                            {"Nitrogen",      0.0},
                            {"Oxygen",        0.0},
                            {"CO2",           0.0},
                            {"Methane",       0.0},
                            {"Argon",         0.0}  // argontha Absolutely Vriliant
                        };
                        recomputeBulkComposition(b);
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
                        int parentIdx = chooseOrbitParent(bodies, worldPos);
                        if (parentIdx >= 0) {
                            Body b;
                            b.mass     = M_EARTH;
                            b.radius   = 6.0;
                            double vol = computeVolumeFromRadius(b.radius);
                            b.density  = b.mass / vol;
                            b.position = worldPos;
                            b.category = BodyCategory::Rocky;
                            b.typeInfo = rockyType;
                            b.albedo   = 0.25;
                            b.name     = "Orbiting body";
                            b.ghost    = true;
                            b.sub = {
                                {"Silicate rock", 0.67},    // 6-7 !!!!!!!!!!! 🎋🎋🌾🌾🪫🪫🥀🥀💔💔😭😭
                                {"Iron",          0.32},
                                {"Water",         0.01},
                                {"Nitrogen",      0.0},
                                {"Oxygen",        0.0},
                                {"CO2",           0.0},
                                {"Methane",       0.0},
                                {"Argon",         0.0}
                            };
                            recomputeBulkComposition(b);

                            Body& parent = bodies[parentIdx];
                            Vec2 r = { worldPos.x - parent.position.x,
                                       worldPos.y - parent.position.y };
                            double dist = length(r);
                            if (dist < SOFTENING_METERS) dist = SOFTENING_METERS;

                            Vec2 rhat = r; normalize(rhat);
                            Vec2 perp{ -rhat.y, rhat.x };
                            double vmag = std::sqrt(G_PHYS * parent.mass / dist); // ts yaabaadoo LMFAO
                            b.velocity = { parent.velocity.x + perp.x * vmag,     // artists who CAN SING vs artists who CANT SING
                                           parent.velocity.y + perp.y * vmag };   // NBA Youngboy CAN SING ✅✅
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
            else if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left) {
                    bodyPanel.dragging = false;
                    controlsPanel.dragging = false;
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
                        b.albedo   = 0.25;
                        b.name     = "Moving body";
                        b.sub = {
                            {"Silicate rock", 0.67},    // SIX SEVEN!!! 🥀🆘🆘🪫🪫🪫💔💔💔💔🧔‍♂️🚡🚡
                            {"Iron",          0.32},
                            {"Water",         0.01},
                            {"Nitrogen",      0.0},
                            {"Oxygen",        0.0},
                            {"CO2",           0.0},
                            {"Methane",       0.0},
                            {"Argon",         0.0}
                        };
                        recomputeBulkComposition(b);
                        b.color    = bodyColor(b);
                        bodies.push_back(b);
                        computeGravity(bodies);
                        updateTemperatures(bodies);
                    }
                }
            }
            // ive lost my shit. i feel like ts is ass but its actually very fine verily i say ergo!
            // "i havent slept properly in weeks" vro tk sh's shkpr 🪫🥀💔🎋🌾
            else if (event.type == sf::Event::MouseMoved) {
                sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
                if (bodyPanel.dragging)
                    bodyPanel.pos = mouseScreenF - bodyPanel.dragOffset;
                if (controlsPanel.dragging)
                    controlsPanel.pos = mouseScreenF - controlsPanel.dragOffset;

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
                    applyOrbitDrag(bodies[orbitEdit.bodyIndex],
                                   bodies[orbitEdit.parentIndex],
                                   mouseWorld);
                }
            }
        } // pollEvent
        // i passed out cold on my desk yesterday (or earlier today depending on your perspective on the matter) and i had a dream about
        // being at an wnba basketball game with santa claus and he opened up his gift-bag and started going, "Ho ho ho! Merry Christmas to all!"
        // while throwing sex toys into the court :sob:. ikyn my brain was making up shit on the spot too, like instead of throwing dildos and
        // vibrators he'd just throw some random ass doohickey like twomad...
        // good thing i wasnt lucid in that cause if i was id be laughing my ass off so hard id probably wake up. i dont get why people get so
        // agitated over shit like this its just goofy as fuck. wbna sucks dick anyhow.
        // ----- update physics -----
        float dtSeconds = clock.restart().asSeconds();
        if (dtSeconds > 0.05f) dtSeconds = 0.05f;
        step(bodies, dtSeconds);

        if (viewBodyIndex >= 0 && viewBodyIndex < (int)bodies.size())
            viewCenterWorld = bodies[viewBodyIndex].position;

        // ----- render -----
        window.clear(sf::Color(10, 10, 20));

        // trails
        if (trailsEnabled) {
            for (const auto& b : bodies) {
                if (b.trail.size() < 2) continue;
                std::vector<sf::Vertex> line;
                line.reserve(b.trail.size());
                for (const auto& p : b.trail)
                    line.emplace_back(worldToScreen(p, WIDTH, HEIGHT),
                                      sf::Color(150, 150, 150));
                window.draw(line.data(), line.size(), sf::LineStrip);
            }
        }

        // AddMoving velocity arrow
        if (isDraggingAddMoving) {
            sf::Vector2f startScreen = worldToScreen(dragStart, WIDTH, HEIGHT);
            sf::Vector2f mouseScreenF = window.mapPixelToCoords(sf::Mouse::getPosition(window));
            sf::Vertex line[] = {
                sf::Vertex(startScreen, sf::Color::White),
                sf::Vertex(mouseScreenF, sf::Color::White)
            };
            window.draw(line, 2, sf::Lines);
        }

        // bodies (with lighting)
        double zoomScale = 1.0e9 / metersPerPixel;
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            const Body& b = bodies[i];

            sf::Vector2f pos = worldToScreen(b.position, WIDTH, HEIGHT);
            float radiusPx = (float)(b.radius * zoomScale);
            if (radiusPx < 2.f)  radiusPx = 2.f;
            if (radiusPx > 60.f) radiusPx = 60.f;

            bool isStarBody = isStar(b);

            // --- Stars: bright core + glow using additive blending ---
            if (isStarBody) {
                sf::Color col = bodyColor(b); // already based on temperature

                // Big outer glow
                float glowOuterR = radiusPx * 2.4f;
                sf::CircleShape glowOuter(glowOuterR);
                glowOuter.setOrigin(glowOuterR, glowOuterR);
                glowOuter.setPosition(pos);
                glowOuter.setFillColor(sf::Color(col.r, col.g, col.b, 40)); // very soft
                window.draw(glowOuter, sf::BlendAdd);

                // Inner glow
                float glowInnerR = radiusPx * 1.6f;
                sf::CircleShape glowInner(glowInnerR);
                glowInner.setOrigin(glowInnerR, glowInnerR);
                glowInner.setPosition(pos);
                glowInner.setFillColor(sf::Color(col.r, col.g, col.b, 90));
                window.draw(glowInner, sf::BlendAdd);

                // Hot core
                float coreR = radiusPx * 0.8f;
                sf::CircleShape core(coreR);
                core.setOrigin(coreR, coreR);
                core.setPosition(pos);

                // core tends toward white to look hotter
                sf::Color coreCol(
                    std::min<int>(255, col.r + 60),
                    std::min<int>(255, col.g + 60),
                    std::min<int>(255, col.b + 60)
                );
                core.setFillColor(coreCol);
                window.draw(core);

                // Selection / fixed / ghost outline on top
                sf::CircleShape outline(radiusPx);
                outline.setOrigin(radiusPx, radiusPx);
                outline.setPosition(pos);
                outline.setFillColor(sf::Color::Transparent);

                if ((int)i == selectedIndex) {
                    outline.setOutlineThickness(2.f);
                    outline.setOutlineColor(sf::Color::Red);
                } else if (b.fixed) {
                    outline.setOutlineThickness(1.5f);
                    outline.setOutlineColor(sf::Color(230, 230, 80));
                } else if (b.ghost) {
                    outline.setOutlineThickness(1.5f);
                    outline.setOutlineColor(sf::Color(0, 150, 255));
                }

                window.draw(outline);
                continue; // skip planet lighting path
            }

            // --- Planets etc with no shader available: flat circle fallback ---
            if (!gHasLightingShader) {
                sf::CircleShape c(radiusPx);
                c.setOrigin(radiusPx, radiusPx);
                c.setPosition(pos);
                c.setFillColor(bodyColor(b));

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
                continue;
            }


            // --- Planets/etc: shaded sprite with shader ---
            // find nearest star for light direction
            int starIdx = findNearestStar(bodies, (int)i);
            Vec2 L{1.0, 0.0};
            float ambient = 0.3f;

            if (starIdx >= 0) {
                Vec2 d{ bodies[starIdx].position.x - b.position.x,
                        bodies[starIdx].position.y - b.position.y }; // planet->star
                double len = length(d);
                if (len > 0.0) {
                    L.x = d.x / len;
                    L.y = d.y / len;
                    ambient = 0.12f;
                }
            }

            // sprite scaled from unit circle texture
            sf::Sprite spr(gPlanetMask);
            spr.setOrigin(128.f, 128.f); // TEX_SIZE / 2
            float scale = radiusPx / 128.f; // radius vs half texture size
            spr.setScale(scale, scale);
            spr.setPosition(pos);

            gPlanetShader.setUniform("texture", gPlanetMask);

            // convert sf::Color -> normalized RGB vec3
            sf::Color c = bodyColor(b);
            sf::Glsl::Vec3 colVec(
                c.r / 255.f,
                c.g / 255.f,
                c.b / 255.f
            );
            gPlanetShader.setUniform("baseColor", colVec);

            gPlanetShader.setUniform("lightDir",
                                    sf::Glsl::Vec2((float)L.x, (float)L.y));
            gPlanetShader.setUniform("ambient", ambient);

            window.draw(spr, &gPlanetShader);

            // outlines (selection / fixed / ghost)
            sf::CircleShape outline(radiusPx);
            outline.setOrigin(radiusPx, radiusPx);
            outline.setPosition(pos);
            outline.setFillColor(sf::Color::Transparent);

            if ((int)i == selectedIndex) {
                outline.setOutlineThickness(2.f);
                outline.setOutlineColor(sf::Color::Red);
            } else if (b.fixed) {
                outline.setOutlineThickness(1.5f);
                outline.setOutlineColor(sf::Color(230, 230, 80));
            } else if (b.ghost) {
                outline.setOutlineThickness(1.5f);
                outline.setOutlineColor(sf::Color(0, 150, 255));
            }

            window.draw(outline);
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
            mText.setPosition(minBtn.getPosition().x + 4.f,
                              minBtn.getPosition().y - 3.f);
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
                l("In body editor: L/R click subcomposition bars");
                l("  to increase/decrease that substance.");
            }
        }

        // Body editor panel + composition UI
        bodyPanelControls.valid = false;
        subCompRects.clear();
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
            minBtn.setPosition(bodyPanel.pos.x + bodyPanel.size.x - 2.f*(btnSize+4.f),
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

            sf::Text mText2;
            mText2.setFont(font);
            mText2.setCharacterSize(12);
            mText2.setFillColor(sf::Color::White);
            mText2.setString(bodyPanel.minimized ? "+" : "-");
            mText2.setPosition(minBtn.getPosition().x + 4.f,
                               minBtn.getPosition().y - 3.f);
            window.draw(mText2);

            if (!bodyPanel.minimized) {
                // --- set up clipping for everything inside the panel ---
                auto winSize = window.getSize();
                unsigned winH = winSize.y;

                // scissor rect is in window coordinates from bottom-left
                GLint sx = static_cast<GLint>(bodyPanel.pos.x);
                GLint sy = static_cast<GLint>(winH - (bodyPanel.pos.y + bodyPanel.size.y));
                GLsizei sw = static_cast<GLsizei>(bodyPanel.size.x);
                GLsizei sh = static_cast<GLsizei>(bodyPanel.size.y);

                glEnable(GL_SCISSOR_TEST);
                glScissor(sx, sy, sw, sh);

                float x = bodyPanel.pos.x + 8.f;
                float y = bodyPanel.pos.y + titleH + 6.0f + bodyPanel.scroll;

                sf::Text line;
                line.setFont(font);
                line.setCharacterSize(13);
                line.setFillColor(sf::Color::White);

                auto drawLine = [&](const std::string& s) {
                    line.setPosition(x, y);
                    line.setString(s);
                    window.draw(line);
                    y += 20.f;
                };

                // --- basic stats ---
                if (renaming)
                    drawLine("Name: " + renameBuffer + "_");
                else
                    drawLine("Name: " + b.name + "  (F2 to rename)");

                drawLine("Type: " + b.typeInfo.name);
                drawLine("Mass: " + formatSci(b.mass) + " kg");
                drawLine("Density: " + formatSci(b.density));
                drawLine("Base radius: " + formatSci(b.radius) + " px");

                if (isStar(b)) {
                    drawLine("Star T: " + formatSci(b.temperatureK) + " K");
                    drawLine("Lum: " + formatSci(b.luminosityW) + " W");
                } else {
                    drawLine("Surface T: " + formatSci(b.temperatureK) + " K");
                    drawLine("Albedo: " + formatSci(b.albedo));
                }

                char buf[128];
                std::snprintf(buf, sizeof(buf),
                            "Comp R/M/I/G/W: %.2f/%.2f/%.2f/%.2f/%.2f",
                            b.bulk.rock, b.bulk.metal, b.bulk.ice,
                            b.bulk.gas, b.bulk.water);
                drawLine(buf);

                // --- primary composition clickable boxes ---
                line.setPosition(x, y);
                line.setString("Primary composition (click):");
                window.draw(line);
                y += 18.f;

                float boxW = 26.f, boxH = 14.f, gap = 6.f;
                float bx = x;
                float by = y;

                auto drawCompBox = [&](const char* label,
                                    const sf::Color& col,
                                    const sf::FloatRect& rect) {
                    sf::RectangleShape r({boxW, boxH});
                    r.setPosition(rect.left, rect.top);
                    r.setFillColor(col);
                    r.setOutlineThickness(1.f);
                    r.setOutlineColor(sf::Color::White);
                    window.draw(r);

                    sf::Text lbl;
                    lbl.setFont(font);
                    lbl.setCharacterSize(11);
                    lbl.setFillColor(sf::Color::White);
                    lbl.setString(label);
                    lbl.setPosition(rect.left, rect.top + boxH + 2.f);
                    window.draw(lbl);
                };

                bodyPanelControls.valid = true;
                bodyPanelControls.compRock  = sf::FloatRect(bx + (boxW+gap)*0, by, boxW, boxH);
                bodyPanelControls.compMetal = sf::FloatRect(bx + (boxW+gap)*1, by, boxW, boxH);
                bodyPanelControls.compIce   = sf::FloatRect(bx + (boxW+gap)*2, by, boxW, boxH);
                bodyPanelControls.compGas   = sf::FloatRect(bx + (boxW+gap)*3, by, boxW, boxH);
                bodyPanelControls.compWater = sf::FloatRect(bx + (boxW+gap)*4, by, boxW, boxH);

                drawCompBox("Rock",  sf::Color(150,120,90), bodyPanelControls.compRock);
                drawCompBox("Metal", sf::Color(180,180,190), bodyPanelControls.compMetal);
                drawCompBox("Ice",   sf::Color(170,220,255), bodyPanelControls.compIce);
                drawCompBox("Gas",   sf::Color(200,170,120), bodyPanelControls.compGas);
                drawCompBox("Water", sf::Color(60,110,170),  bodyPanelControls.compWater);

                y += boxH + 26.f;

                // --- detailed subcomposition bars ---
                line.setPosition(x, y);
                line.setString("Detailed composition (L/R click bars):");
                window.draw(line);
                y += 18.f;

                subCompRects.clear();
                const float barW = bodyPanel.size.x - 24.f;
                const float barH = 10.f;

                for (std::size_t i = 0; i < b.sub.size() && i < 7; ++i) {
                    const auto& sc = b.sub[i];
                    float frac = (float)std::clamp(sc.fraction, 0.0, 1.0);

                    sf::RectangleShape bg({barW, barH});
                    bg.setPosition(x, y);
                    bg.setFillColor(sf::Color(40, 40, 60));
                    window.draw(bg);

                    sf::RectangleShape fg({barW * frac, barH});
                    fg.setPosition(x, y);
                    fg.setFillColor(sf::Color(100, 180, 220));
                    window.draw(fg);

                    sf::Text lbl;
                    lbl.setFont(font);
                    lbl.setCharacterSize(11);
                    lbl.setFillColor(sf::Color::White);
                    std::snprintf(buf, sizeof(buf), "%s: %.2f", sc.name.c_str(), sc.fraction);
                    lbl.setString(buf);
                    lbl.setPosition(x, y + barH + 1.f);
                    window.draw(lbl);

                    subCompRects.emplace_back(x, y, barW, barH);
                    y += barH + 18.f;
                }

                std::string flags = "Fixed: ";
                flags += (b.fixed ? "Yes" : "No");
                flags += "   Volume lock: ";
                flags += (lockVolume ? "On" : "Off");
                drawLine(flags);

                glDisable(GL_SCISSOR_TEST);
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