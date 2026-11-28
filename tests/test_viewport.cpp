// ============================================================================
// Tests: Viewport (viewport.h): zoom and pan
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <vector>

#include "viewport.h"
#include "test_framework.h"
#include "test_support.h"

// --- Zoom and pan (viewport) ---
TEST(viewport_zooms_about_the_cursor) {
    Viewport v(1000.0f, 1000.0f);
    CHECK(v.zoom() == 1.0f && v.width() == 1000.0f);
    // Zooming in about (300, 700): that point stays at the same place on screen (30 % across, 70 % down)
    v.zoomAt(300.0f, 700.0f, 4.0f);
    CHECK(v.zoom() == 4.0f && v.width() == 250.0f);
    CHECK(std::fabs((300.0f - v.left()) / v.width() - 0.3f) < 1e-5f);
    CHECK(std::fabs((700.0f - v.top()) / v.height() - 0.7f) < 1e-5f);
    // Can't zoom out past the whole map, or in past maxZoom()
    v.zoomAt(500.0f, 500.0f, 0.01f);
    CHECK(v.zoom() == 1.0f && v.centreX() == 500.0f && v.centreY() == 500.0f);
    v.zoomAt(500.0f, 500.0f, 1000.0f);
    CHECK(v.zoom() == v.maxZoom());
}

TEST(viewport_stays_inside_the_map) {
    Viewport v(1000.0f, 1000.0f);
    v.zoomAt(10.0f, 10.0f, 5.0f);                    // Zoom near a corner
    CHECK(v.left() >= 0.0f && v.top() >= 0.0f);
    v.pan(-5000.0f, 5000.0f);                         // Pan far off the map
    CHECK(v.left() == 0.0f && std::fabs(v.top() + v.height() - 1000.0f) < 1e-3f);
    v.reset();
    CHECK(v.zoom() == 1.0f && v.left() == 0.0f);
}

TEST(viewport_frames_the_fire) {
    std::vector<std::uint8_t> mask(N, 0);
    int x0, y0, x1, y1;
    CHECK(!maskBounds(mask, GRID_WIDTH, GRID_HEIGHT, x0, y0, x1, y1));
    for (int y = 240; y <= 260; ++y)
        for (int x = 230; x <= 250; ++x) mask[at(x, y)] = 1;
    CHECK(maskBounds(mask, GRID_WIDTH, GRID_HEIGHT, x0, y0, x1, y1));
    CHECK(x0 == 230 && x1 == 250 && y0 == 240 && y1 == 260);

    // Frame a 21-cell fire with 10 cells of margin: 41 cells across, centred on the fire
    Viewport v(static_cast<float>(GRID_WIDTH), static_cast<float>(GRID_HEIGHT));   // 1 world unit per cell
    v.frame(230.0f, 240.0f, 251.0f, 261.0f, 10.0f);
    // The view keeps the map's aspect ratio, so on a non-square map only the tighter side is exactly 41 cells;
    // the other is wider, and the whole 41 x 41 box is always in view
    CHECK(std::fabs(std::min(v.width(), v.height()) - 41.0f) < 1e-3f);
    CHECK(v.width() >= 41.0f - 1e-3f && v.height() >= 41.0f - 1e-3f);
    CHECK(std::fabs(v.centreX() - 240.5f) < 1e-3f && std::fabs(v.centreY() - 250.5f) < 1e-3f);
}
