#pragma once
#include <cstdint>
#include <vector>

// ============================================================================
// Viewport: which part of the map is on screen (zoom and pan)
//
// Pure maths in map ("world") units, so it is part of wf_core and testable without SFML. The renderer turns
// it into an sf::View. World units are window pixels at zoom 1 (cells x CELL_SIZE).
// ============================================================================

constexpr float VIEW_MAX_ZOOM = 50.0f;   // At most 50x: 10 cells across a 500-cell map
constexpr float VIEW_ZOOM_STEP = 1.25f;  // One mouse-wheel notch or +/- key press

class Viewport {
public:
    Viewport(float mapWidth, float mapHeight);

    float centreX() const { return cx_; }
    float centreY() const { return cy_; }
    float zoom() const { return zoom_; }                 // 1 = whole map visible
    float width() const { return mapW_ / zoom_; }        // Visible size, world units
    float height() const { return mapH_ / zoom_; }
    float left() const { return cx_ - width() / 2; }
    float top() const { return cy_ - height() / 2; }

    // Zooms by `factor` (> 1 = in) keeping the world point (wx, wy) under the cursor where it is on screen
    void zoomAt(float wx, float wy, float factor);
    // Moves the view by (dx, dy) world units
    void pan(float dx, float dy);
    // Shows the box [minX, maxX] x [minY, maxY] (world units) with `margin` around it, as large as possible
    void frame(float minX, float minY, float maxX, float maxY, float margin);
    // Back to the whole map
    void reset();

private:
    void clamp();   // Keeps the zoom in range and the view inside the map

    float mapW_, mapH_;
    float cx_, cy_;
    float zoom_ = 1.0f;
};

// Bounding box (in cells, inclusive) of the non-zero entries of a width x height mask. False if all zero.
bool maskBounds(const std::vector<std::uint8_t>& mask, int width, int height,
                int& minX, int& minY, int& maxX, int& maxY);
