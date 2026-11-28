// ============================================================================
// Tests: Wind arrow overlay (wind_overlay.h)
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "fire_model.h"
#include "viewport.h"
#include "wind_overlay.h"
#include "wind_solver.h"
#include "test_framework.h"
#include "test_support.h"

// ----------------------------------------------------------------------------
// Wind arrow overlay (wind_overlay.h)
// ----------------------------------------------------------------------------
namespace {
    // A wind field that is the same everywhere: speedKmh, blowing FROM fromDeg
    WindField uniformWind(float speedKmh, float fromDeg) {
        WindField f;
        f.grid = G;
        f.toward.assign(N, windTowardGrid(fromDeg));
        f.speedKmh.assign(N, speedKmh);
        return f;
    }

    // World units per screen pixel for a window as wide as the map (CELL_SIZE pixels per cell at zoom 1)
    float perPixelAt(const Viewport& v) { return v.width() / static_cast<float>(GRID_WIDTH * CELL_SIZE); }
}

// Arrows stay about WIND_ARROW_SPACING_PX apart on screen at every zoom, cover the view, and stay put in map
// space while panning (the lattice is fixed to the map)
TEST(wind_arrows_keep_their_screen_spacing_and_stay_put_when_panning) {
    const WindField f = uniformWind(20.0f, 270.0f);
    for (const float zoom : { 1.0f, 3.0f, 8.0f, 40.0f }) {
        Viewport v(static_cast<float>(GRID_WIDTH * CELL_SIZE), static_cast<float>(GRID_HEIGHT * CELL_SIZE));
        v.zoomAt(v.centreX(), v.centreY(), zoom);
        const float perPixel = perPixelAt(v);
        const float stepPx = windArrowSpacing(G, perPixel) / perPixel;   // Lattice spacing on screen
        const std::vector<WindArrow> arrows = buildWindArrows(f, v, perPixel, 20.0f, 20.0f);
        // Every arrow's centre lies in the view or within one spacing of it
        const float step = windArrowSpacing(G, perPixel);
        bool inView = true;
        for (const WindArrow& a : arrows) {
            const float mx = 0.5f * (a.x0 + a.x1), my = 0.5f * (a.y0 + a.y1);
            inView = inView && mx >= v.left() - step && mx <= v.left() + v.width() + step
                            && my >= v.top() - step && my <= v.top() + v.height() + step;
        }
        // About (view / spacing)^2 arrows: the view is covered
        const float across = v.width() / step;
        std::printf("  (zoom %.0f: spacing %.0f px, %zu arrows)\n", zoom, stepPx, arrows.size());
        CHECK(stepPx >= WIND_ARROW_SPACING_PX && stepPx < 2.0f * WIND_ARROW_SPACING_PX + 1e-3f);
        CHECK(inView);
        CHECK(arrows.size() >= static_cast<std::size_t>(across * across * 0.8f));
    }

    // Pan by a third of a spacing: arrows present in both views are at exactly the same map positions
    Viewport v(static_cast<float>(GRID_WIDTH * CELL_SIZE), static_cast<float>(GRID_HEIGHT * CELL_SIZE));
    v.zoomAt(v.centreX(), v.centreY(), 6.0f);
    const float perPixel = perPixelAt(v);
    const std::vector<WindArrow> before = buildWindArrows(f, v, perPixel, 20.0f, 20.0f);
    v.pan(windArrowSpacing(G, perPixel) / 3.0f, 0.0f);
    const std::vector<WindArrow> after = buildWindArrows(f, v, perPixel, 20.0f, 20.0f);
    int shared = 0, moved = 0;
    for (const WindArrow& a : after) {
        bool found = false;
        for (const WindArrow& b : before) found = found || (a.x0 == b.x0 && a.y0 == b.y0 && a.x1 == b.x1 && a.y1 == b.y1);
        (found ? shared : moved) += 1;
    }
    std::printf("  (after a pan: %d arrows unchanged, %d new at the edge)\n", shared, moved);
    CHECK(shared > 0);
    CHECK(moved <= static_cast<int>(after.size()) / 5);   // Only a new column at the edge
}

// Each arrow points where the wind blows, with its length in proportion to the speed and its colour from the
// speed against the open wind: white = same, blue = slower, magenta = faster
TEST(wind_arrows_show_direction_speed_and_colour) {
    Viewport v(static_cast<float>(GRID_WIDTH * CELL_SIZE), static_cast<float>(GRID_HEIGHT * CELL_SIZE));
    const float perPixel = perPixelAt(v);
    const float step = windArrowSpacing(G, perPixel);

    // Open wind from the west everywhere: arrows point east (+x), white, 85 % of the spacing long
    const std::vector<WindArrow> east = buildWindArrows(uniformWind(20.0f, 270.0f), v, perPixel, 20.0f, 20.0f);
    CHECK(!east.empty());
    const WindArrow& a = east.front();
    CHECK(a.x1 > a.x0 && std::fabs(a.y1 - a.y0) < 1e-3f);
    CHECK(std::fabs((a.x1 - a.x0) - 0.85f * step) < 1e-3f);
    CHECK(a.hx0 < a.x1 && a.hx1 < a.x1);                     // Head strokes trail back from the tip...
    CHECK((a.hy0 - a.y1) * (a.hy1 - a.y1) < 0.0f);           // ...one either side of the shaft
    CHECK(a.r == 255 && a.g == 255 && a.b == 255);

    // Wind from the south (blows north, -y), half the fastest speed: half the length
    const std::vector<WindArrow> north = buildWindArrows(uniformWind(10.0f, 180.0f), v, perPixel, 20.0f, 20.0f);
    CHECK(north.front().y1 < north.front().y0);
    CHECK(std::fabs((north.front().y0 - north.front().y1) - 0.425f * step) < 1e-3f);

    // Colours: 30 % or more slower = fully blue, 30 % or more faster = fully magenta, in between = partly coloured
    std::uint8_t r, g, b;
    windArrowColour(14.0f, 20.0f, r, g, b);
    CHECK(r == 60 && g == 140 && b == 255);
    windArrowColour(26.0f, 20.0f, r, g, b);
    CHECK(r == 255 && g == 60 && b == 200);
    windArrowColour(23.0f, 20.0f, r, g, b);
    CHECK(r == 255 && g > 60 && g < 255 && b > 200 && b < 255);

    // No field (USE_TERRAIN_WIND off, or nothing solved): no arrows
    CHECK(buildWindArrows(WindField(), v, perPixel, 20.0f, 20.0f).empty());
}

// On the terrain wind over a hill, the arrow at the crest is the longest and magenta (faster than the open
// wind), and the one at the upwind foot is shorter and blue (slower): the overlay shows what the solver found
TEST(wind_arrows_show_the_terrain_wind_over_a_hill) {
    WindSolverSettings s;
    s.columnSizeM = 600.0f;
    const WindField f = solveTerrainWind(gridWithCells(WIND_TEST_CELL_M), gaussianHill(), 20.0f, 270.0f, s);
    const float maxKmh = *std::max_element(f.speedKmh.begin(), f.speedKmh.end());
    Viewport v(static_cast<float>(GRID_WIDTH * CELL_SIZE), static_cast<float>(GRID_HEIGHT * CELL_SIZE));
    const float perPixel = perPixelAt(v);
    const std::vector<WindArrow> arrows = buildWindArrows(f, v, perPixel, 20.0f, maxKmh);
    // The arrows nearest the crest and nearest the upwind foot (two standard deviations west of it)
    const float cs = static_cast<float>(CELL_SIZE);
    const auto nearest = [&](float x, float y) {
        const WindArrow* best = nullptr;
        float bestD = 1e30f;
        for (const WindArrow& a : arrows) {
            const float dx = 0.5f * (a.x0 + a.x1) - x, dy = 0.5f * (a.y0 + a.y1) - y;
            if (dx * dx + dy * dy < bestD) { bestD = dx * dx + dy * dy; best = &a; }
        }
        return *best;
    };
    const WindArrow crest = nearest(CX * cs, CY * cs);
    const WindArrow foot = nearest((CX - 2.0f * 2000.0f / WIND_TEST_CELL_M) * cs, CY * cs);
    const auto length = [](const WindArrow& a) { return std::hypot(a.x1 - a.x0, a.y1 - a.y0); };
    std::printf("  (crest %.1f km/h colour %d,%d,%d; upwind foot %.1f km/h colour %d,%d,%d)\n",
                crest.speedKmh, crest.r, crest.g, crest.b, foot.speedKmh, foot.r, foot.g, foot.b);
    CHECK(crest.speedKmh > 20.0f && foot.speedKmh < 20.0f);
    CHECK(length(crest) > length(foot));
    CHECK(crest.g < 255 && crest.r == 255);   // Magenta side
    CHECK(foot.r < 255 && foot.b == 255);     // Blue side
}
