#include "renderer.h"
#include <algorithm>
#include <cmath>

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
        "2D Wildfire Simulator (Rothermel, LANDFIRE) - wheel: zoom, right-drag: pan, F: frame fire"),
    terrainShade_(terrainShade),
    pixels_(static_cast<size_t>(GRID_WIDTH)* GRID_HEIGHT * 4, 255),       // RGBA primary grid buffer
    overlayPixels_(static_cast<size_t>(GRID_WIDTH)* GRID_HEIGHT * 4, 0),   // RGBA mission overlay buffer
    marker_(3.0f), sensor_(SENSOR_RADIUS* CELL_SIZE), body_(5.0f),
    view_(static_cast<float>(GRID_WIDTH * CELL_SIZE), static_cast<float>(GRID_HEIGHT * CELL_SIZE)) {
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
    applyView();
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
// - Unburned: Grayscale elevation shading (white = highest elevation); road cells tinted amber.
// - Burning: orange-red surface fire, orange torching, yellow active crown fire.
// - Burned: Dark red-brown char color.
void Renderer::paintBase(int idx, const std::vector<Cell>& grid) {
    switch (grid[idx].state) {
    case CellState::Unburned: {
        const std::uint8_t g = terrainShade_[idx];
        if (showRoads && grid[idx].roadClass != 0) {
            // Roads: amber, blended half-and-half with the elevation shade so the terrain still shows through
            setColor(idx, static_cast<std::uint8_t>((g + 230) / 2), static_cast<std::uint8_t>((g + 160) / 2),
                     static_cast<std::uint8_t>(g / 2));
        }
        else {
            setColor(idx, g, g, g);
        }
        break;
    }
    case CellState::Burning:
        // Colour by fire type so crown fire runs stand out (as FlamMap's crown fire activity map)
        if (grid[idx].fireType == 2)      setColor(idx, 255, 230, 60);   // Active crown fire: yellow
        else if (grid[idx].fireType == 1) setColor(idx, 255, 150, 0);    // Torching (passive): orange
        else                              setColor(idx, 255, 69, 0);     // Surface fire: orange-red
        break;
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
    // World units per screen pixel: markers are scaled by this so they keep their on-screen size when zoomed
    const float perPixel = view_.width() / static_cast<float>(window_.getSize().x);

    // Cell grid lines, once cells are big enough on screen to tell apart (fading in as you zoom)
    const float cellPx = cs / perPixel;
    if (cellPx >= 8.0f) {
        const sf::Color lineColour(0, 0, 0, static_cast<sf::Uint8>(std::min(80.0f, 20.0f + (cellPx - 8.0f) * 4.0f)));
        const int x0 = std::max(0, static_cast<int>(view_.left() / cs));
        const int x1 = std::min(GRID_WIDTH, static_cast<int>(std::ceil((view_.left() + view_.width()) / cs)));
        const int y0 = std::max(0, static_cast<int>(view_.top() / cs));
        const int y1 = std::min(GRID_HEIGHT, static_cast<int>(std::ceil((view_.top() + view_.height()) / cs)));
        sf::VertexArray lines(sf::Lines);
        for (int x = x0; x <= x1; ++x) {   // Vertical lines on cell edges
            lines.append(sf::Vertex(sf::Vector2f(x * cs, y0 * cs), lineColour));
            lines.append(sf::Vertex(sf::Vector2f(x * cs, y1 * cs), lineColour));
        }
        for (int y = y0; y <= y1; ++y) {   // Horizontal lines
            lines.append(sf::Vertex(sf::Vector2f(x0 * cs, y * cs), lineColour));
            lines.append(sf::Vertex(sf::Vector2f(x1 * cs, y * cs), lineColour));
        }
        window_.draw(lines);
    }
    marker_.setScale(perPixel, perPixel);
    body_.setScale(perPixel, perPixel);

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
// ============================================================================
// Zoom and Pan
// ============================================================================

// Copies the viewport into the window. Everything drawn afterwards (grid, overlays, drones) uses it, and
// mapPixelToCoords (mouse clicks) follows it too.
void Renderer::applyView() {
    sf::View v(sf::FloatRect(view_.left(), view_.top(), view_.width(), view_.height()));
    window_.setView(v);
}

// ----------------------------------------------------------------------------
// Mouse and keyboard controls for the view
// ----------------------------------------------------------------------------
bool Renderer::handleViewEvent(const sf::Event& event) {
    const float perPixel = view_.width() / static_cast<float>(window_.getSize().x);

    switch (event.type) {
    case sf::Event::MouseWheelScrolled: {
        // Zoom about the point under the cursor: one notch = VIEW_ZOOM_STEP
        if (event.mouseWheelScroll.wheel != sf::Mouse::VerticalWheel) return false;
        const sf::Vector2f at = window_.mapPixelToCoords(sf::Vector2i(event.mouseWheelScroll.x, event.mouseWheelScroll.y));
        view_.zoomAt(at.x, at.y, std::pow(VIEW_ZOOM_STEP, event.mouseWheelScroll.delta));
        applyView();
        return true;
    }
    case sf::Event::MouseButtonPressed:
        // Right or middle button starts a drag (left button is left free for igniting cells)
        if (event.mouseButton.button != sf::Mouse::Right && event.mouseButton.button != sf::Mouse::Middle) return false;
        dragging_ = true;
        lastMouse_ = sf::Vector2i(event.mouseButton.x, event.mouseButton.y);
        return true;
    case sf::Event::MouseButtonReleased:
        if (event.mouseButton.button != sf::Mouse::Right && event.mouseButton.button != sf::Mouse::Middle) return false;
        dragging_ = false;
        return true;
    case sf::Event::MouseMoved: {
        if (!dragging_) return false;
        // The map follows the mouse: move the view the opposite way to the cursor
        const sf::Vector2i now(event.mouseMove.x, event.mouseMove.y);
        view_.pan((lastMouse_.x - now.x) * perPixel, (lastMouse_.y - now.y) * perPixel);
        lastMouse_ = now;
        applyView();
        return true;
    }
    case sf::Event::KeyPressed: {
        const float step = 0.1f * view_.width();   // Arrow keys move a tenth of the visible width
        switch (event.key.code) {
        case sf::Keyboard::Add: case sf::Keyboard::Equal:
            view_.zoomAt(view_.centreX(), view_.centreY(), VIEW_ZOOM_STEP); break;
        case sf::Keyboard::Subtract: case sf::Keyboard::Hyphen:
            view_.zoomAt(view_.centreX(), view_.centreY(), 1.0f / VIEW_ZOOM_STEP); break;
        case sf::Keyboard::Left:  view_.pan(-step, 0); break;
        case sf::Keyboard::Right: view_.pan(step, 0); break;
        case sf::Keyboard::Up:    view_.pan(0, -step); break;
        case sf::Keyboard::Down:  view_.pan(0, step); break;
        case sf::Keyboard::Num0: case sf::Keyboard::Numpad0: case sf::Keyboard::Home:
            view_.reset(); break;
        default:
            return false;
        }
        applyView();
        return true;
    }
    default:
        return false;
    }
}

// ----------------------------------------------------------------------------
// Frame the fire: fit the bounding box of the mask, with 10 cells of margin
// ----------------------------------------------------------------------------
void Renderer::frameMask(const std::vector<std::uint8_t>& mask) {
    int x0, y0, x1, y1;
    if (!maskBounds(mask, GRID_WIDTH, GRID_HEIGHT, x0, y0, x1, y1)) return;   // No fire yet
    const float cs = static_cast<float>(CELL_SIZE);
    view_.frame(x0 * cs, y0 * cs, (x1 + 1) * cs, (y1 + 1) * cs, 10.0f * cs);
    applyView();
}
