#include "renderer.h"

namespace {
    // Color scheme for different drone mission path types
    const sf::Color kSweepColour(60, 140, 255); // Soft blue for sector sweep drones
    const sf::Color kAcoColour(255, 80, 255);   // Magenta for ACO outer perimeter tour drones
}

// ============================================================================
// Real-time Visual Renderer (SFML Graphics Engine)
//
// Visualizes:
// 1. Terrain elevation shading and dynamic wildfire combustion states.
// 2. Dynamic containment perimeter rings (cyan).
// 3. Voronoi sector boundaries, standoff safety zones, and outer search limits.
// 4. Multi-drone positions, sensor range radii, and active flight paths.
// ============================================================================
Renderer::Renderer(const std::vector<std::uint8_t>& terrainShade)
    : window_(sf::VideoMode(GRID_WIDTH* CELL_SIZE, GRID_HEIGHT* CELL_SIZE),
        "2D Wildfire Simulator (GDAL GIS Engine) - P: perimeter, V: drone mission overlay"),
    terrainShade_(terrainShade),
    pixels_(static_cast<size_t>(GRID_WIDTH)* GRID_HEIGHT * 4, 255),       // RGBA primary grid buffer
    overlayPixels_(static_cast<size_t>(GRID_WIDTH)* GRID_HEIGHT * 4, 0),   // RGBA mission overlay buffer
    marker_(3.0f), sensor_(SENSOR_RADIUS* CELL_SIZE), body_(5.0f) {
    window_.setFramerateLimit(30); // Cap frame rate to 30 FPS

    // Create primary terrain & fire texture/sprite mapped to grid dimensions
    texture_.create(GRID_WIDTH, GRID_HEIGHT);
    texture_.update(pixels_.data());
    sprite_.setTexture(texture_);
    sprite_.setScale(static_cast<float>(CELL_SIZE), static_cast<float>(CELL_SIZE));

    // Create mission overlay texture/sprite (Voronoi borders, boundaries)
    overlayTex_.create(GRID_WIDTH, GRID_HEIGHT);
    overlaySprite_.setTexture(overlayTex_);
    overlaySprite_.setScale(static_cast<float>(CELL_SIZE), static_cast<float>(CELL_SIZE));

    // Configure shape graphics for POI markers, drone sensors, and drone bodies
    marker_.setOrigin(3.0f, 3.0f);
    sensor_.setOrigin(SENSOR_RADIUS * CELL_SIZE, SENSOR_RADIUS * CELL_SIZE);
    sensor_.setFillColor(sf::Color::Transparent);
    sensor_.setOutlineThickness(1.0f);
    body_.setOrigin(5.0f, 5.0f);
    body_.setOutlineColor(sf::Color::White);
    body_.setOutlineThickness(1.0f);
}

// Helper function to set raw RGBA pixel values in the base texture pixel array
void Renderer::setColor(int idx, std::uint8_t r, std::uint8_t g, std::uint8_t b) {
    pixels_[idx * 4 + 0] = r;
    pixels_[idx * 4 + 1] = g;
    pixels_[idx * 4 + 2] = b;
}

// ----------------------------------------------------------------------------
// Cell Color Mapping
// ----------------------------------------------------------------------------
// Assigns pixel RGB colors according to current cell combustion state:
// - Unburned: Grayscale elevation shading (white = highest elevation).
// - Burning: Bright orange/red fire color.
// - Burned: Dark red-brown char color.
void Renderer::paintBase(int idx, const std::vector<Cell>& grid) {
    switch (grid[idx].state) {
    case CellState::Unburned: {
        const std::uint8_t g = terrainShade_[idx];
        setColor(idx, g, g, g);
        break;
    }
    case CellState::Burning: setColor(idx, 255, 69, 0); break;   // Orange-Red
    case CellState::Burned:  setColor(idx, 110, 35, 25); break;  // Charcoal / Red-brown
    }
}

// Repaints the entire base grid texture array
void Renderer::paintAll(const std::vector<Cell>& grid) {
    for (int i = 0; i < GRID_WIDTH * GRID_HEIGHT; ++i) paintBase(i, grid);
}

// Repaints only dirty/modified grid cells for performance optimization
void Renderer::paintCells(const std::vector<Cell>& grid, const std::vector<int>& cells) {
    for (int idx : cells) paintBase(idx, grid);
}

// Sets or clears cyan containment perimeter ring highlights
void Renderer::setRing(const std::vector<int>& ring, const std::vector<Cell>& grid) {
    for (int idx : ringCells_) paintBase(idx, grid); // Restore previous ring cells to base terrain color
    ringCells_.clear();

    if (!showPerimeter) return;

    for (int idx : ring) setColor(idx, 0, 255, 255); // Highlight perimeter cells in Cyan
    ringCells_ = ring;
}

// ----------------------------------------------------------------------------
// Mission Visual Overlay Generator
// ----------------------------------------------------------------------------
// Renders:
// 1. Voronoi sector boundary lines (White).
// 2. Standoff safety zone edge (Red).
// 3. Outer search zone limit (Yellow).
// 4. Closed flight paths for all drones.
void Renderer::setMissionVisuals(const MissionPlan& plan, const std::vector<Drone>& drones,
    const DistanceTransform& liveEdt) {
    std::fill(overlayPixels_.begin(), overlayPixels_.end(), 0); // Reset overlay pixel buffer

    auto labelAt = [&](int x, int y) {
        return (x < 0 || x >= GRID_WIDTH || y < 0 || y >= GRID_HEIGHT) ? -1 : plan.label[y * GRID_WIDTH + x];
        };
    auto holeAt = [&](int x, int y) {
        return x >= 0 && x < GRID_WIDTH && y >= 0 && y < GRID_HEIGHT
            && liveEdt.squaredDistance(y * GRID_WIDTH + x) < NO_FLY_SQ;
        };

    // Draw lines along cell borders in overlay texture
    for (int y = 0; y < GRID_HEIGHT; ++y) {
        for (int x = 0; x < GRID_WIDTH; ++x) {
            const int l = labelAt(x, y);
            if (l < 0) continue;

            const int r = labelAt(x + 1, y), dn = labelAt(x, y + 1),
                lf = labelAt(x - 1, y), up = labelAt(x, y - 1);

            sf::Uint8* px = &overlayPixels_[(y * GRID_WIDTH + x) * 4];

            if (r < 0 || dn < 0 || lf < 0 || up < 0) { // Boundary edge of active search zone
                if (holeAt(x + 1, y) || holeAt(x, y + 1) || holeAt(x - 1, y) || holeAt(x, y - 1)) {
                    // Standoff edge adjacent to fire (Red)
                    px[0] = 255; px[1] = 60; px[2] = 60; px[3] = 220;
                }
                else {
                    // Outer search limit (Yellow)
                    px[0] = 255; px[1] = 220; px[2] = 0; px[3] = 200;
                }
            }
            else if (r != l || dn != l) { // Interior Voronoi region border (White)
                px[0] = 255; px[1] = 255; px[2] = 255; px[3] = 110;
            }
        }
    }
    overlayTex_.update(overlayPixels_.data());

    // Generate line-strip path polylines for each drone's target route
    paths_.assign(drones.size(), sf::VertexArray(sf::LineStrip));
    for (size_t i = 0; i < drones.size(); ++i) {
        const Drone& d = drones[i];
        sf::Color c = d.isACO ? kAcoColour : kSweepColour;
        c.a = 120; // Semi-transparent path line

        const float cs = static_cast<float>(CELL_SIZE);
        for (const Vec2& p : d.route) {
            paths_[i].append(sf::Vertex(sf::Vector2f(p.x * cs, p.y * cs), c));
        }

        // Loop path back to start point to close the circuit
        if (d.route.size() > 1) {
            paths_[i].append(sf::Vertex(sf::Vector2f(d.route.front().x * cs, d.route.front().y * cs), c));
        }
    }
}

// ----------------------------------------------------------------------------
// Render Frame Execution
// ----------------------------------------------------------------------------
// Displays base terrain, fire state updates, mission overlays, POIs, and dynamic 
// drone entities onto the SFML display window.
void Renderer::draw(const MissionPlan& plan, const std::vector<Drone>& drones) {
    texture_.update(pixels_.data());
    window_.clear(sf::Color::Black);
    window_.draw(sprite_);

    const float cs = static_cast<float>(CELL_SIZE);

    // Draw mission overlay details if enabled ('V' key toggle)
    if (showMission) {
        window_.draw(overlaySprite_);
        for (const sf::VertexArray& path : paths_) window_.draw(path);

        // Draw region POI markers (White = Inner sector centroid, Yellow = Outer ACO tour POI)
        for (const Region& rg : plan.regions) {
            if (rg.innerCount > 0) {
                marker_.setFillColor(sf::Color::White);
                marker_.setPosition(rg.innerCentroid.x * cs, rg.innerCentroid.y * cs);
                window_.draw(marker_);
            }
            if (rg.outerCount > 0) {
                marker_.setFillColor(sf::Color(255, 220, 0));
                marker_.setPosition(rg.poi.x * cs, rg.poi.y * cs);
                window_.draw(marker_);
            }
        }
    }

    // Render drone sensor coverage bubbles and physical drone body markers
    for (const Drone& d : drones) {
        const sf::Color c = d.isACO ? kAcoColour : kSweepColour;

        // Render sensor radius circle outline
        sensor_.setOutlineColor(sf::Color(c.r, c.g, c.b, 160));
        sensor_.setPosition(d.pos.x * cs, d.pos.y * cs);
        window_.draw(sensor_);

        // Render drone body marker
        body_.setFillColor(c);
        body_.setPosition(d.pos.x * cs, d.pos.y * cs);
        window_.draw(body_);
    }

    window_.display();
}