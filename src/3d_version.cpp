#include <SFML/Graphics.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>

// ----------------------
// Basic 3D vector
// ----------------------
struct Vec3 {
    double x, y, z;

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s)     const { return {x * s,   y * s,   z * s};   }

    Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vec3& operator*=(double s)      { x *= s;   y *= s;   z *= s;   return *this; }
    Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
};

double length(const Vec3& v) {
    return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}
void normalize(Vec3& v) {
    double len = length(v);
    if (len > 0.0) {
        v.x /= len; v.y /= len; v.z /= len;
    }
}
double dot(const Vec3& a, const Vec3& b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

bool gridEnabled = true;

// ----------------------
// Body
// ----------------------
struct Body {
    double mass;        // kg (or arbitrary units)
    Vec3   position;    // m (or arbitrary units)
    Vec3   velocity;    // m/s
    Vec3   acceleration;

    double radius;      // visual radius in "world units"
    sf::Color color;

    std::string name;
    bool fixed = false;

    std::vector<Vec3> trail;
};

// ----------------------
// Physics constants / globals
// ----------------------
const double G_PHYS = 6.67430e-11;  // real G (you can fudge it if you like)
double timeScale = 86400.0;         // 1 day of sim per real second

bool trailsEnabled = true;
int  maxTrailPoints = 300;
bool paused = false;

// ----------------------
// Camera
// ----------------------
struct Camera {
    Vec3   target{0.0, 0.0, 0.0};  // point we orbit around
    double distance = 5.0e10;      // distance from target (world units)
    double yaw   = 0.8;            // radians (rotation around Y)
    double pitch = 0.4;            // radians (up/down)
};

Camera camera;

// Compute camera position from yaw/pitch/distance
Vec3 getCameraPosition() {
    double cy = std::cos(camera.yaw);
    double sy = std::sin(camera.yaw);
    double cp = std::cos(camera.pitch);
    double sp = std::sin(camera.pitch);

    Vec3 offset {
        camera.distance * cp * cy,
        camera.distance * sp,
        camera.distance * cp * sy
    };
    return camera.target + offset;
}

void getCameraBasis(Vec3& forward, Vec3& right, Vec3& up) {
    Vec3 camPos = getCameraPosition();

    forward = { camera.target.x - camPos.x,
                camera.target.y - camPos.y,
                camera.target.z - camPos.z };
    normalize(forward);

    Vec3 worldUp{0.0, 1.0, 0.0};
    right = cross(worldUp, forward);
    if (length(right) < 1e-8) {
        right = {1.0, 0.0, 0.0};
    } else {
        normalize(right);
    }
    up = cross(forward, right);
}

// Forward declaration so drawGridPlane can use it
sf::Vector2f projectToScreen(const Vec3& world,
                             unsigned width, unsigned height,
                             bool& onScreen);

void drawGridPlane(sf::RenderWindow& window,
                   unsigned width, unsigned height)
{
    if (!gridEnabled) return;

    const double extent = 3.0e11;
    const double step   = 5.0e10;

    sf::Color minorColor(30, 30, 60);
    sf::Color majorColor(60, 60, 100);
    
    void drawGridPlane(sf::RenderWindow& window, unsigned width, unsigned height);
    {
        if (!gridEnabled) return;
    
        // Total half-size of grid in world units (centered at origin on z=0)
        const double extent = 3.0e11;   // ~2 AU radius if your planet is at 1.5e11
        const double step   = 5.0e10;   // spacing between minor grid lines
    
        sf::Color minorColor(30, 30, 60);
        sf::Color majorColor(60, 60, 100);
    
        auto drawGridLine = [&](const Vec3& A, const Vec3& B, sf::Color col) {
            bool onA = false, onB = false;
            sf::Vector2f a2 = projectToScreen(A, width, height, onA);
            sf::Vector2f b2 = projectToScreen(B, width, height, onB);
            if (!onA && !onB) return;
    
            sf::Vertex line[] = {
                sf::Vertex(a2, col),
                sf::Vertex(b2, col)
            };
            window.draw(line, 2, sf::Lines);
        };
    
        // Lines parallel to X (vary Y)
        for (double y = -extent; y <= extent + 1.0; y += step) {
            bool major = (std::fmod(std::fabs(y), step * 5.0) < 1.0e-3);
            sf::Color col = major ? majorColor : minorColor;
    
            Vec3 p1{-extent, y, 0.0};
            Vec3 p2{+extent, y, 0.0};
            drawGridLine(p1, p2, col);
        }
    
        // Lines parallel to Y (vary X)
        for (double x = -extent; x <= extent + 1.0; x += step) {
            bool major = (std::fmod(std::fabs(x), step * 5.0) < 1.0e-3);
            sf::Color col = major ? majorColor : minorColor;
    
            Vec3 p1{x, -extent, 0.0};
            Vec3 p2{x, +extent, 0.0};
            drawGridLine(p1, p2, col);
        }
    }
}


// Project a 3D point onto the screen
sf::Vector2f projectToScreen(const Vec3& world,
                             unsigned width, unsigned height,
                             bool& onScreen)
{
    Vec3 camPos = getCameraPosition();

    // Build camera basis
    Vec3 forward = {camera.target.x - camPos.x,
                    camera.target.y - camPos.y,
                    camera.target.z - camPos.z};
    normalize(forward);

    Vec3 worldUp{0.0, 1.0, 0.0};
    Vec3 right = cross(worldUp, forward);
    if (length(right) < 1e-8) {
        right = {1.0, 0.0, 0.0};
    } else {
        normalize(right);
    }
    Vec3 up = cross(forward, right);

    // Convert point to camera space
    Vec3 rel {
        world.x - camPos.x,
        world.y - camPos.y,
        world.z - camPos.z
    };

    double xCam = dot(rel, right);
    double yCam = dot(rel, up);
    double zCam = dot(rel, forward); // "depth"

    if (zCam <= 0.0) {
        onScreen = false;
        return {0.f, 0.f};
    }

    // Perspective projection
    double fov = 60.0 * 3.141592653589793 / 180.0;
    double f = 0.5 * height / std::tan(fov * 0.5);

    double xNDC = (xCam * f) / zCam;
    double yNDC = (yCam * f) / zCam;

    float cx = width * 0.5f;
    float cy = height * 0.5f;

    sf::Vector2f screenPos(
        cx + static_cast<float>(xNDC),
        cy - static_cast<float>(yNDC)
    );
    onScreen = true;
    return screenPos;
}

// ----------------------
// Gravity & integration
// ----------------------
void computeGravity(std::vector<Body>& bodies) {
    std::size_t n = bodies.size();
    for (auto& b : bodies) {
        b.acceleration = {0.0, 0.0, 0.0};
    }

    const double softening2 = 1.0e18; // soften small distances

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            Vec3 r = bodies[j].position - bodies[i].position;
            double dist2 = r.x*r.x + r.y*r.y + r.z*r.z + softening2;
            double dist  = std::sqrt(dist2);
            double invDist3 = 1.0 / (dist2 * dist);

            Vec3 acc_i = r * (G_PHYS * bodies[j].mass * invDist3);
            Vec3 acc_j = r * (-G_PHYS * bodies[i].mass * invDist3);

            if (!bodies[i].fixed) bodies[i].acceleration += acc_i;
            if (!bodies[j].fixed) bodies[j].acceleration += acc_j;
        }
    }
}

void step(std::vector<Body>& bodies, double dtRealSeconds) {
    if (paused) return;

    // dtRealSeconds is the real frame time;
    // timeScale is sim seconds per real second
    double dtSim = dtRealSeconds * timeScale;

    const double MAX_SUBSTEP = 3600.0; // 1 simulated hour max per substep
    double absDt = std::abs(dtSim);
    int nSub = (absDt <= MAX_SUBSTEP) ? 1
                                      : (int)std::ceil(absDt / MAX_SUBSTEP);
    double h = dtSim / nSub;

    bool recordTrails = trailsEnabled;
    if (!recordTrails) {
        for (auto& b : bodies)
            b.trail.clear();
    }

    for (int s = 0; s < nSub; ++s) {
        // v_{n+1} = v_n + a_n * h
        for (auto& b : bodies) {
            if (!b.fixed) {
                b.velocity += b.acceleration * h;
            }
        }
        // x_{n+1} = x_n + v_{n+1} * h
        for (auto& b : bodies) {
            if (!b.fixed) {
                b.position += b.velocity * h;
            }
        }

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

// ----------------------
// Utility
// ----------------------
int findNearestBody3D(const std::vector<Body>& bodies,
                      const Vec3& worldPos,
                      double maxDist)
{
    int bestIndex = -1;
    double best2 = maxDist * maxDist;
    for (std::size_t i = 0; i < bodies.size(); ++i) {
        Vec3 d = bodies[i].position - worldPos;
        double d2 = d.x*d.x + d.y*d.y + d.z*d.z;
        if (d2 < best2) {
            best2 = d2;
            bestIndex = (int)i;
        }
    }
    return bestIndex;
}

std::string formatSci(double v, int prec = 2) {
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.*e", prec, v);
    return std::string(buf);
}

// ----------------------
// Main
// ----------------------
int main() {
    const unsigned WIDTH = 1280;
    const unsigned HEIGHT = 720;

    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT),
                            "PlanetSim3D.exe");
    window.setFramerateLimit(120);

    sf::Font font;
    bool fontLoaded = font.loadFromFile("resources/arial.ttf");

    // ---------------------------------
    // Build an initial simple system
    // ---------------------------------
    std::vector<Body> bodies;

    // Star at origin
    Body star;
    star.mass     = 1.98847e30;  // Sun
    star.position = {0.0, 0.0, 0.0};
    star.velocity = {0.0, 0.0, 0.0};
    star.radius   = 7.0e9;       // visual radius (arbitrary)
    star.color    = sf::Color(255, 240, 180);
    star.name     = "Star";
    star.fixed    = true;
    bodies.push_back(star);

    // One planet in XY plane
    Body planet;
    planet.mass     = 5.97219e24;   // Earth
    planet.position = {1.5e11, 0.0, 0.0}; // ~1 AU along +X
    // circular orbit speed v = sqrt(GM/r)
    double rOrbit = length(planet.position - star.position);
    double vOrbit = std::sqrt(G_PHYS * star.mass / rOrbit);
    planet.velocity = {0.0, vOrbit, 0.0}; // around +Y for a simple orbit
    planet.radius   = 3.0e9;              // visual radius
    planet.color    = sf::Color(80, 140, 220);
    planet.name     = "Planet";
    planet.fixed    = false;
    bodies.push_back(planet);

    computeGravity(bodies);

    int selectedIndex = -1;

    // Camera initial setup (looking at origin)
    camera.target   = {0.0, 0.0, 0.0};
    camera.distance = 4.0e11; // far enough to see orbit
    camera.yaw      = 0.7;
    camera.pitch    = 0.35;

    // For adding bodies: project mouse into a reference plane (z = 0)
    enum class InteractionMode { AddStill, AddMoving };
    InteractionMode mode = InteractionMode::AddStill;
    bool draggingAddMoving = false;
    Vec3 dragStartWorld{0,0,0};

    sf::Clock clock;

    while (window.isOpen()) {
        // -----------------
        // Event handling
        // -----------------
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            else if (event.type == sf::Event::KeyPressed) {
                if (event.key.code == sf::Keyboard::Space) {
                    paused = !paused;
                }
                else if (event.key.code == sf::Keyboard::Num1) {
                    mode = InteractionMode::AddStill;
                }
                else if (event.key.code == sf::Keyboard::Num2) {
                    mode = InteractionMode::AddMoving;
                }
                else if (event.key.code == sf::Keyboard::Up) {
                    camera.pitch += 0.05;
                    if (camera.pitch > 1.4) camera.pitch = 1.4;
                }
                else if (event.key.code == sf::Keyboard::Down) {
                    camera.pitch -= 0.05;
                    if (camera.pitch < -1.4) camera.pitch = -1.4;
                }
                else if (event.key.code == sf::Keyboard::Left) {
                    camera.yaw -= 0.05;
                }
                else if (event.key.code == sf::Keyboard::Right) {
                    camera.yaw += 0.05;
                }
                // Zoom with R / F
                else if (event.key.code == sf::Keyboard::F) {
                    camera.distance *= 0.8;
                    if (camera.distance < 1.0e10) camera.distance = 1.0e10;
                }
                else if (event.key.code == sf::Keyboard::R) {
                    camera.distance *= 1.25;
                    if (camera.distance > 1.0e13) camera.distance = 1.0e13;
                }
                else if (event.key.code == sf::Keyboard::LBracket) {
                    timeScale *= 0.5;
                    if (timeScale < 3600.0) timeScale = 3600.0;
                }
                else if (event.key.code == sf::Keyboard::RBracket) {
                    timeScale *= 2.0;
                    if (timeScale > 1.0e7) timeScale = 1.0e7;
                }
                else if (event.key.code == sf::Keyboard::G) {
                    gridEnabled = !gridEnabled;
                }
                else if (event.key.code == sf::Keyboard::T) {
                    trailsEnabled = !trailsEnabled;
                }
                else if (event.key.code == sf::Keyboard::Delete) {
                    if (selectedIndex > 0 &&
                        selectedIndex < (int)bodies.size()) {
                        bodies.erase(bodies.begin() + selectedIndex);
                        selectedIndex = -1;
                        computeGravity(bodies);
                    }
                }
            }
            else if (event.type == sf::Event::MouseButtonPressed) {
                sf::Vector2i mousePix = sf::Mouse::getPosition(window);

                // Raycast into z=0 plane
                Vec3 camPos = getCameraPosition();

                Vec3 forward = {camera.target.x - camPos.x,
                                camera.target.y - camPos.y,
                                camera.target.z - camPos.z};
                normalize(forward);
                Vec3 worldUp{0.0, 1.0, 0.0};
                Vec3 right = cross(worldUp, forward);
                if (length(right) < 1e-8) right = {1.0,0.0,0.0};
                else normalize(right);
                Vec3 up = cross(forward, right);

                float cx = WIDTH * 0.5f;
                float cy = HEIGHT * 0.5f;
                float sx = (float)mousePix.x - cx;
                float sy = cy - (float)mousePix.y;

                double fov = 60.0 * 3.141592653589793 / 180.0;
                double f = 0.5 * HEIGHT / std::tan(fov * 0.5);

                Vec3 dirCam{ sx, sy, f };
                normalize(dirCam);

                Vec3 dirWorld{
                    right.x * dirCam.x + up.x * dirCam.y + forward.x * dirCam.z,
                    right.y * dirCam.x + up.y * dirCam.y + forward.y * dirCam.z,
                    right.z * dirCam.x + up.z * dirCam.y + forward.z * dirCam.z
                };
                normalize(dirWorld);

                double t = 0.0;
                if (std::fabs(dirWorld.z) > 1e-9) {
                    t = -camPos.z / dirWorld.z;
                }
                if (t < 0.0) t = 0.0;

                Vec3 worldPos {
                    camPos.x + dirWorld.x * t,
                    camPos.y + dirWorld.y * t,
                    camPos.z + dirWorld.z * t
                };

                if (event.mouseButton.button == sf::Mouse::Left) {
                    if (mode == InteractionMode::AddStill) {
                        Body b;
                        b.mass     = 5.0e24;
                        b.position = worldPos;
                        b.velocity = {0.0, 0.0, 0.0};
                        b.radius   = 2.0e9;
                        b.color    = sf::Color(160, 120, 90);
                        b.name     = "Body";
                        bodies.push_back(b);
                        computeGravity(bodies);
                    } else if (mode == InteractionMode::AddMoving) {
                        draggingAddMoving = true;
                        dragStartWorld = worldPos;
                    }
                }
                else if (event.mouseButton.button == sf::Mouse::Right) {
                    int idx = findNearestBody3D(bodies, worldPos, 5.0e10);
                    selectedIndex = idx;
                    if (idx >= 0) {
                        camera.target = bodies[idx].position;
                    }
                }
            }
            else if (event.type == sf::Event::MouseButtonReleased) {
                if (event.mouseButton.button == sf::Mouse::Left &&
                    draggingAddMoving) {

                    draggingAddMoving = false;

                    sf::Vector2i mousePix = sf::Mouse::getPosition(window);

                    Vec3 camPos = getCameraPosition();
                    Vec3 forward = {camera.target.x - camPos.x,
                                    camera.target.y - camPos.y,
                                    camera.target.z - camPos.z};
                    normalize(forward);
                    Vec3 worldUp{0.0, 1.0, 0.0};
                    Vec3 right = cross(worldUp, forward);
                    if (length(right) < 1e-8) right = {1.0,0.0,0.0};
                    else normalize(right);
                    Vec3 up = cross(forward, right);

                    float cx = WIDTH * 0.5f;
                    float cy = HEIGHT * 0.5f;
                    float sx = (float)mousePix.x - cx;
                    float sy = cy - (float)mousePix.y;

                    double fov = 60.0 * 3.141592653589793 / 180.0;
                    double ff = 0.5 * HEIGHT / std::tan(fov * 0.5);

                    Vec3 dirCam{ sx, sy, ff };
                    normalize(dirCam);
                    Vec3 dirWorld{
                        right.x * dirCam.x + up.x * dirCam.y + forward.x * dirCam.z,
                        right.y * dirCam.x + up.y * dirCam.y + forward.y * dirCam.z,
                        right.z * dirCam.x + up.z * dirCam.y + forward.z * dirCam.z
                    };
                    normalize(dirWorld);

                    double t = 0.0;
                    if (std::fabs(dirWorld.z) > 1e-9) {
                        t = -camPos.z / dirWorld.z;
                    }
                    if (t < 0.0) t = 0.0;

                    Vec3 worldEnd{
                        camPos.x + dirWorld.x * t,
                        camPos.y + dirWorld.y * t,
                        camPos.z + dirWorld.z * t
                    };

                    Vec3 delta = worldEnd - dragStartWorld;

                    Body b;
                    b.mass     = 5.0e24;
                    b.position = dragStartWorld;
                    b.velocity = delta * 1.0e-4;
                    b.radius   = 2.0e9;
                    b.color    = sf::Color(120, 200, 150);
                    b.name     = "Moving body";

                    bodies.push_back(b);
                    computeGravity(bodies);
                }
            }
        } // end pollEvent

        // -----------------
        // Update physics
        // -----------------
        float dtReal = clock.restart().asSeconds();
        if (dtReal > 0.05f) dtReal = 0.05f;
        step(bodies, dtReal);

        // -----------------
        // Camera movement (WASDQE)
        // -----------------
        Vec3 forwardBasis, rightBasis, upBasis;
        getCameraBasis(forwardBasis, rightBasis, upBasis);

        Vec3 moveDir{0.0, 0.0, 0.0};

        if (sf::Keyboard::isKeyPressed(sf::Keyboard::W)) {
            moveDir += forwardBasis;       // forward
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::S)) {
            moveDir -= forwardBasis;       // backward
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::A)) {
            moveDir -= rightBasis;         // left
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::D)) {
            moveDir += rightBasis;         // right
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Q)) {
            moveDir += upBasis;            // up
        }
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::E)) {
            moveDir -= upBasis;            // down
        }

        double lenMove = length(moveDir);
        if (lenMove > 0.0) {
            moveDir *= (1.0 / lenMove);

            double baseSpeed = camera.distance * 0.4;
            if (sf::Keyboard::isKeyPressed(sf::Keyboard::LShift) ||
                sf::Keyboard::isKeyPressed(sf::Keyboard::RShift)) {
                baseSpeed *= 3.0;
            }

            double stepDist = baseSpeed * dtReal;
            Vec3 delta = moveDir * stepDist;

            camera.target += delta;
        }

        // If following a body, keep camera.target synced
        if (selectedIndex >= 0 && selectedIndex < (int)bodies.size()) {
            camera.target = bodies[selectedIndex].position;
        }

        // -----------------
        // Render
        // -----------------
        window.clear(sf::Color(5, 5, 15));

        drawGridPlane(window, WIDTH, HEIGHT);

        // Trails
        if (trailsEnabled) {
            for (const auto& b : bodies) {
                if (b.trail.size() < 2) continue;
                std::vector<sf::Vertex> verts;
                verts.reserve(b.trail.size());
                for (const auto& p3 : b.trail) {
                    bool onScreen = false;
                    sf::Vector2f p2 = projectToScreen(p3, WIDTH, HEIGHT, onScreen);
                    if (!onScreen) continue;
                    verts.emplace_back(p2, sf::Color(150,150,150));
                }
                if (verts.size() >= 2)
                    window.draw(verts.data(), verts.size(), sf::LineStrip);
            }
        }

        // Bodies
        Vec3 camPos = getCameraPosition();
        for (std::size_t i = 0; i < bodies.size(); ++i) {
            const Body& b = bodies[i];

            bool onScreen = false;
            sf::Vector2f pos2 = projectToScreen(b.position, WIDTH, HEIGHT, onScreen);
            if (!onScreen) continue;

            Vec3 diff = b.position - camPos;
            double dist = length(diff);
            if (dist <= 0.0) dist = 1.0;

            double baseDist = 1.0e11;
            double scale = baseDist / dist;
            float radiusPx = (float)(b.radius * scale / 5.0e9);
            if (radiusPx < 2.0f)  radiusPx = 2.0f;
            if (radiusPx > 40.0f) radiusPx = 40.0f;

            sf::CircleShape circ(radiusPx);
            circ.setOrigin(radiusPx, radiusPx);
            circ.setPosition(pos2);
            circ.setFillColor(b.color);

            if ((int)i == selectedIndex) {
                circ.setOutlineThickness(2.0f);
                circ.setOutlineColor(sf::Color::Red);
            } else if (b.fixed) {
                circ.setOutlineThickness(1.5f);
                circ.setOutlineColor(sf::Color(230, 230, 80));
            }

            window.draw(circ);
        }

        // HUD
        if (fontLoaded) {
            sf::Text text;
            text.setFont(font);
            text.setCharacterSize(14);
            text.setFillColor(sf::Color::White);

            std::string modeStr = (mode == InteractionMode::AddStill)
                                  ? "AddStill(1)"
                                  : "AddMoving(2)";
            std::string sel = "None";
            if (selectedIndex >= 0 && selectedIndex < (int)bodies.size())
                sel = bodies[selectedIndex].name;

            char buf[128];
            std::snprintf(buf, sizeof(buf),
                          "Mode: %s   Selected: %s   timeScale = %.2e s/s   %s",
                          modeStr.c_str(), sel.c_str(), timeScale,
                          paused ? "[PAUSED]" : "");
            text.setString(buf);
            text.setPosition(10.f, 10.f);
            window.draw(text);

            text.setCharacterSize(13);
            std::snprintf(buf, sizeof(buf),
                          "Camera: dist=%.2e  yaw=%.2f  pitch=%.2f",
                          camera.distance, camera.yaw, camera.pitch);
            text.setString(buf);
            text.setPosition(10.f, HEIGHT - 24.f);
            window.draw(text);

            std::string controls =
                "Space: pause   [ / ]: timeScale\n"
                "Arrows: yaw/pitch   R/F: zoom in/out\n"
                "W/S: forward/back   A/D: left/right\n"
                "Q/E: up/down   Shift: faster\n"
                "1: Add still   2: Add moving\n"
                "LMB: add body / drag velocity\n"
                "RMB: select + follow";
            text.setString(controls);
            text.setPosition(WIDTH - 280.f, HEIGHT - 80.f);
            window.draw(text);
        }

        window.display();
    }
    return 0;
}