#include "viewport.h"
#include <algorithm>

// ============================================================================
// Viewport (see viewport.h)
// ============================================================================

Viewport::Viewport(float mapWidth, float mapHeight)
    : mapW_(mapWidth), mapH_(mapHeight), cx_(mapWidth / 2), cy_(mapHeight / 2) {}

// ----------------------------------------------------------------------------
// Zoom about a point
// ----------------------------------------------------------------------------
// The point keeps the same fraction across the view before and after, so it stays under the mouse cursor.
void Viewport::zoomAt(float wx, float wy, float factor) {
    const float fx = (wx - left()) / width();     // Where the point sits across the view, 0..1
    const float fy = (wy - top()) / height();
    zoom_ = std::clamp(zoom_ * factor, 1.0f, VIEW_MAX_ZOOM);
    cx_ = wx - fx * width() + width() / 2;
    cy_ = wy - fy * height() + height() / 2;
    clamp();
}

void Viewport::pan(float dx, float dy) {
    cx_ += dx;
    cy_ += dy;
    clamp();
}

// ----------------------------------------------------------------------------
// Fit a box (e.g. the fire) in the view
// ----------------------------------------------------------------------------
void Viewport::frame(float minX, float minY, float maxX, float maxY, float margin) {
    const float bw = (maxX - minX) + 2 * margin, bh = (maxY - minY) + 2 * margin;
    zoom_ = std::min(mapW_ / std::max(bw, 1e-3f), mapH_ / std::max(bh, 1e-3f));   // Largest zoom that fits it
    cx_ = (minX + maxX) / 2;
    cy_ = (minY + maxY) / 2;
    clamp();
}

void Viewport::reset() {
    zoom_ = 1.0f;
    cx_ = mapW_ / 2;
    cy_ = mapH_ / 2;
}

// ----------------------------------------------------------------------------
// Keep the zoom between 1x and VIEW_MAX_ZOOM, and never show anything outside the map
// ----------------------------------------------------------------------------
void Viewport::clamp() {
    zoom_ = std::clamp(zoom_, 1.0f, VIEW_MAX_ZOOM);
    cx_ = std::clamp(cx_, width() / 2, mapW_ - width() / 2);
    cy_ = std::clamp(cy_, height() / 2, mapH_ - height() / 2);
}

// ----------------------------------------------------------------------------
// Bounding box of a mask (used to frame the fire)
// ----------------------------------------------------------------------------
bool maskBounds(const std::vector<std::uint8_t>& mask, int width, int height,
                int& minX, int& minY, int& maxX, int& maxY) {
    minX = width; minY = height; maxX = -1; maxY = -1;
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (!mask[static_cast<size_t>(y) * width + x]) continue;
            minX = std::min(minX, x); maxX = std::max(maxX, x);
            minY = std::min(minY, y); maxY = std::max(maxY, y);
        }
    }
    return maxX >= 0;
}
