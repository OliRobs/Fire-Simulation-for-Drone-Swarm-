#include "wind_overlay.h"
#include <algorithm>
#include <cmath>

// ============================================================================
// Wind arrow overlay (what each arrow means is in wind_overlay.h)
// ============================================================================

namespace {
    constexpr float FULL_LENGTH = 0.85f;     // Fastest wind's arrow, as a fraction of the lattice spacing
    constexpr float MIN_LENGTH = 0.15f;      // Shortest arrow drawn, same units, so calm cells still show a stub
    constexpr float HEAD_FRACTION = 0.35f;   // Arrowhead stroke length / arrow length
    constexpr float HEAD_MIN = 0.12f;        // ...but at least this fraction of the spacing, so stubs keep a head
    constexpr float HEAD_ANGLE_RAD = 0.49f;  // ~28 deg between each head stroke and the shaft
    // 30 % slower / faster than the open wind = fully blue / magenta. Terrain changes the 20-ft wind by roughly
    // -35 % to +25 % on the San Gabriel and fallback maps, so this spreads the colours over that range.
    constexpr float FULL_COLOUR_AT = 0.3f;

    // Linear blend of two 0-255 channels
    std::uint8_t blend(int from, int to, float t) {
        return static_cast<std::uint8_t>(std::lround(from + (to - from) * t));
    }
}

// ----------------------------------------------------------------------------
// Colour from the speed against the open wind: white at the open wind's speed, blending to blue (slower) or
// magenta (faster), fully coloured at FULL_COLOUR_AT (30 %) slower or faster
// ----------------------------------------------------------------------------
void windArrowColour(float speedKmh, float openKmh, std::uint8_t& r, std::uint8_t& g, std::uint8_t& b) {
    r = g = b = 255;
    if (openKmh <= 0.0f) return;   // No open wind to compare with: plain white
    const float ratio = speedKmh / openKmh;
    if (ratio < 1.0f) {            // Slower: towards blue (60, 140, 255)
        const float t = std::min(1.0f, (1.0f - ratio) / FULL_COLOUR_AT);
        r = blend(255, 60, t);
        g = blend(255, 140, t);
    }
    else if (ratio > 1.0f) {       // Faster: towards magenta (255, 60, 200)
        const float t = std::min(1.0f, (ratio - 1.0f) / FULL_COLOUR_AT);
        g = blend(255, 60, t);
        b = blend(255, 200, t);
    }
}

// ----------------------------------------------------------------------------
// Lattice spacing: whole powers of two of cells, so the lattice points of a finer spacing include those of
// every coarser one, and the arrows only re-arrange when the zoom crosses a factor of two
// ----------------------------------------------------------------------------
float windArrowSpacing(float perPixel, float spacingPx) {
    const float cellsWanted = spacingPx * perPixel / static_cast<float>(CELL_SIZE);   // Cells per spacingPx on screen
    int cells = 1;
    while (cells < cellsWanted && cells < std::max(GRID_WIDTH, GRID_HEIGHT)) cells *= 2;
    return static_cast<float>(cells * CELL_SIZE);
}

// ----------------------------------------------------------------------------
// One arrow per lattice point inside the view, centred on the point
// ----------------------------------------------------------------------------
std::vector<WindArrow> buildWindArrows(const WindField& field, const Viewport& view, float perPixel,
                                       float openKmh, float maxKmh, float spacingPx) {
    std::vector<WindArrow> arrows;
    const std::size_t cellCount = static_cast<std::size_t>(GRID_WIDTH) * GRID_HEIGHT;
    if (field.toward.size() != cellCount || field.speedKmh.size() != cellCount || maxKmh <= 0.0f) return arrows;

    const float step = windArrowSpacing(perPixel, spacingPx);            // World units between arrows
    const float cs = static_cast<float>(CELL_SIZE);
    // Lattice points sit at (k + 0.5) x step; take those inside the view plus up to one step beyond each edge,
    // so an arrow whose centre is just off screen but whose shaft reaches in is still drawn.
    // (k + 0.5) step >= left - step  <=>  k >= left / step - 1.5;  (k + 0.5) step <= right + step  <=>  k <= right / step + 0.5
    const int kx0 = std::max(0, static_cast<int>(std::ceil(view.left() / step - 1.5f)));
    const int ky0 = std::max(0, static_cast<int>(std::ceil(view.top() / step - 1.5f)));
    const int kx1 = static_cast<int>(std::floor((view.left() + view.width()) / step + 0.5f));
    const int ky1 = static_cast<int>(std::floor((view.top() + view.height()) / step + 0.5f));
    const float cosH = std::cos(HEAD_ANGLE_RAD), sinH = std::sin(HEAD_ANGLE_RAD);

    for (int ky = ky0; ky <= ky1; ++ky) {
        const float py = (ky + 0.5f) * step;
        const int cy = static_cast<int>(py / cs);   // Cell under the lattice point
        if (cy < 0 || cy >= GRID_HEIGHT) continue;
        for (int kx = kx0; kx <= kx1; ++kx) {
            const float px = (kx + 0.5f) * step;
            const int cx = static_cast<int>(px / cs);
            if (cx < 0 || cx >= GRID_WIDTH) continue;
            const std::size_t idx = static_cast<std::size_t>(cy) * GRID_WIDTH + cx;
            const Vec2 dir = field.toward[idx];   // Unit vector, grid frame (x east, y south) = world frame
            const float speed = field.speedKmh[idx];

            WindArrow a;
            a.speedKmh = speed;
            // Length in proportion to the speed; centred on the lattice point
            const float len = step * std::max(MIN_LENGTH, FULL_LENGTH * std::min(1.0f, speed / maxKmh));
            a.x0 = px - 0.5f * len * dir.x;
            a.y0 = py - 0.5f * len * dir.y;
            a.x1 = px + 0.5f * len * dir.x;
            a.y1 = py + 0.5f * len * dir.y;
            // Arrowhead: two strokes back from the tip, the reversed direction turned by +-HEAD_ANGLE_RAD
            const float head = std::min(len, std::max(HEAD_FRACTION * len, HEAD_MIN * step));
            const float bx = -dir.x, by = -dir.y;
            a.hx0 = a.x1 + head * (bx * cosH - by * sinH);
            a.hy0 = a.y1 + head * (bx * sinH + by * cosH);
            a.hx1 = a.x1 + head * (bx * cosH + by * sinH);
            a.hy1 = a.y1 + head * (-bx * sinH + by * cosH);
            windArrowColour(speed, openKmh, a.r, a.g, a.b);
            arrows.push_back(a);
        }
    }
    return arrows;
}
