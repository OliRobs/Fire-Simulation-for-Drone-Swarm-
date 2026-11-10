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
void computeTopography(std::vector<Cell>& grid, float cellSizeX, float cellSizeY, bool haveSlopeLayer) {
    constexpr float RAD_TO_DEG = 57.2957795f;

    // Elevation at (x, y), clamped at the map edge so edge cells use a one-sided difference
    auto z = [&](int x, int y) {
        x = std::clamp(x, 0, GRID_WIDTH - 1);
        y = std::clamp(y, 0, GRID_HEIGHT - 1);
        return grid[y * GRID_WIDTH + x].elevation;
    };

    for (int y = 0; y < GRID_HEIGHT; ++y) {
        for (int x = 0; x < GRID_WIDTH; ++x) {
            const float a = z(x - 1, y - 1), b = z(x, y - 1), c = z(x + 1, y - 1);
            const float d = z(x - 1, y),                     f = z(x + 1, y);
            const float g = z(x - 1, y + 1), h = z(x, y + 1), i = z(x + 1, y + 1);

            const float dzdx = ((c + 2 * f + i) - (a + 2 * d + g)) / (8.0f * cellSizeX);
            const float dzdy = ((g + 2 * h + i) - (a + 2 * b + c)) / (8.0f * cellSizeY);
            const float tanSlope = std::sqrt(dzdx * dzdx + dzdy * dzdy);

            Cell& cell = grid[y * GRID_WIDTH + x];
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
