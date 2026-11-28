#include "topography.h"
#include <algorithm>
#include <cmath>

// ============================================================================
// Slope and Aspect from the DEM
//
// LANDFIRE SlpD gives slope steepness, but the Rothermel vector addition also needs which way is uphill.
// That comes from the elevation gradient, estimated with Horn's (1981) 3x3 weighted differences, the method
// GIS packages (and LANDFIRE's own slope/aspect products) use:
//
//     a b c      dz/dx = ((c + 2f + i) - (a + 2d + g)) / (8 * cellSizeX)    (x = east)
//     d e f      dz/dy = ((g + 2h + i) - (a + 2b + c)) / (8 * cellSizeY)    (y = south, i.e. down the rows)
//     g h i
//
// The gradient (dz/dx, dz/dy) points uphill; its length is tan(slope).
// ============================================================================
void computeTopography(const GridSpec& spec, std::vector<Cell>& grid, bool haveSlopeLayer) {
    constexpr float RAD_TO_DEG = 57.2957795f;
    requireGridSize(grid, spec, "computeTopography cells");
    const int W = spec.width, H = spec.height;
    const float cellSizeX = spec.cellSizeX, cellSizeY = spec.cellSizeY;

    // Elevation at (x, y), clamped at the map edge so edge cells use a one-sided difference
    auto z = [&](int x, int y) {
        x = std::clamp(x, 0, W - 1);
        y = std::clamp(y, 0, H - 1);
        return grid[y * W + x].elevation;
    };

    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            const float a = z(x - 1, y - 1), b = z(x, y - 1), c = z(x + 1, y - 1);
            const float d = z(x - 1, y),                     f = z(x + 1, y);
            const float g = z(x - 1, y + 1), h = z(x, y + 1), i = z(x + 1, y + 1);

            const float dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8.0f * cellSizeX);
            const float dzdy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8.0f * cellSizeY);
            const float tanSlope = std::sqrt(dzdx * dzdx + dzdy * dzdy);

            Cell& cell = grid[y * W + x];
            if (tanSlope > 1e-6f) {
                cell.upslopeX = dzdx / tanSlope;
                cell.upslopeY = dzdy / tanSlope;
            }
            else {
                cell.upslopeX = cell.upslopeY = 0.0f;   // Flat: no uphill direction
            }
            if (!haveSlopeLayer) cell.slopeDeg = std::atan(tanSlope) * RAD_TO_DEG;
        }
    }
}

// ============================================================================
// Simulation-cell terrain from a finer DEM (USGS 3DEP 1 m)
//
// Each fine pixel gets its own Horn gradient. Per simulation cell:
//   elevation = mean of the pixel elevations
//   slope     = mean of the pixel slopes (degrees). Averaging the slopes rather than the elevations keeps
//               the real steepness of rough ground, which matters because Rothermel's slope factor grows
//               with tan^2(slope) (RMRS-GTR-371, phi_s = 5.275 beta^-0.3 tan^2(slope)).
//   upslope   = direction of the mean gradient vector (sum of the pixel gradients), i.e. the overall uphill
//               direction across the cell. On a sharp ridge the two sides partly cancel, as they should.
// ============================================================================
void aggregateFineTopography(const std::vector<float>& fine, int cellsX, int cellsY, int factor,
                             float finePixelX, float finePixelY, std::vector<CellTopography>& out) {
    constexpr double RAD_TO_DEG = 57.29577951308232;
    const int fineW = cellsX * factor + 2;   // Buffer width including the halo

    // Running sums per cell (double: a 10 m cell sums 100 pixels, a 1 km strip far more)
    struct Sums { double elev = 0, slope = 0, gx = 0, gy = 0; int n = 0; };
    std::vector<Sums> sums(static_cast<size_t>(cellsX) * cellsY);

    for (int py = 0; py < cellsY * factor; ++py) {
        const int by = py + 1;                        // Row in the buffer (skipping the top halo row)
        const int cellRow = py / factor;
        for (int px = 0; px < cellsX * factor; ++px) {
            const int bx = px + 1;
            const float e = fine[static_cast<size_t>(by) * fineW + bx];
            if (std::isnan(e)) continue;              // No data: this pixel adds nothing

            // Neighbour elevation, falling back to the centre where the neighbour has no data
            auto z = [&](int dx, int dy) {
                const float v = fine[static_cast<size_t>(by + dy) * fineW + (bx + dx)];
                return std::isnan(v) ? e : v;
            };
            const float a = z(-1, -1), b = z(0, -1), c = z(1, -1);
            const float d = z(-1, 0),                f = z(1, 0);
            const float g = z(-1, 1),  h = z(0, 1),  i = z(1, 1);

            // Horn (1981) 3x3 gradient, as in computeTopography (x east, y south)
            const double dzdx = ((c + 2.0 * f + i) - (a + 2.0 * d + g)) / (8.0 * finePixelX);
            const double dzdy = ((g + 2.0 * h + i) - (a + 2.0 * b + c)) / (8.0 * finePixelY);

            Sums& s = sums[static_cast<size_t>(cellRow) * cellsX + px / factor];
            s.elev += e;
            s.slope += std::atan(std::sqrt(dzdx * dzdx + dzdy * dzdy)) * RAD_TO_DEG;
            s.gx += dzdx;
            s.gy += dzdy;
            ++s.n;
        }
    }

    // Turn the sums into cell values
    out.assign(sums.size(), CellTopography{});
    for (size_t k = 0; k < sums.size(); ++k) {
        const Sums& s = sums[k];
        CellTopography& t = out[k];
        t.validPixels = s.n;
        if (s.n == 0) continue;                       // Not covered by the DEM
        t.elevation = static_cast<float>(s.elev / s.n);
        t.slopeDeg = static_cast<float>(s.slope / s.n);
        const double gx = s.gx / s.n, gy = s.gy / s.n, len = std::sqrt(gx * gx + gy * gy);
        if (len > 1e-6) {                             // Same flat-ground threshold as computeTopography
            t.upslopeX = static_cast<float>(gx / len);
            t.upslopeY = static_cast<float>(gy / len);
        }
    }
}
