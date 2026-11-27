#pragma once
#include <vector>
#include "config.h"
#include "types.h"

// ============================================================================
// Terrain-adjusted wind: a mass-conserving ("COM") wind solver, after WindNinja
//
// The fire model needs a 20-ft wind in every cell. Weather gives one open wind (WIND_SPEED_20FT_KMH from
// WIND_FROM_DEG); real wind speeds up over ridges, is channelled along valleys and turns round hills. This
// module adjusts the one wind over the terrain the way WindNinja's default solver does (Forthofer et al. 2014,
// Int. J. Wildland Fire 23: 969-981; github.com/firelab/windninja, src/ninja/ninja.cpp):
//
//  1. Mesh. The fire grid is grouped into square columns of about columnSizeM (10 x 10 fire cells at the
//     default 100 m). Each column becomes one stack of `layers` mesh nodes, from the column's mean ground
//     elevation up to a flat top topAboveSummitM above the highest summit, thin near the ground and growing by
//     layerGrowth per layer (terrain-following). WindNinja puts its top much lower (see WIND_TOP_ABOVE_SUMMIT_M
//     in config.h for why this one is higher). The 3D wind is solved on this coarse mesh, not on the 10 m fire cells:
//     at 10 m it would be ~24 million unknowns, and ridge/valley effects are hundreds of metres across.
//  2. First guess (u0, v0, w0). The open wind is spread up each column with a neutral log profile,
//     U(z) = U20 ln((z - d) / z0) / ln((6.1 m + h - d) / z0), where h is the vegetation height (canopy, or the
//     fuel bed depth where there are no trees), z0 = 0.13 h and d = 0.64 h. Those two ratios are the ones the
//     Albini & Baughman wind adjustment factor in rothermel.h is built on, so the two parts agree. Vertical
//     wind starts at 0. This field ignores the terrain, so it piles air into hillsides (it has divergence).
//  3. Adjustment. Find the field (u, v, w) closest to the first guess that conserves mass (div u = 0):
//         minimise  integral of  (u - u0)^2 + (v - v0)^2 + (alpha_V/alpha_H)^2 (w - w0)^2   with div u = 0.
//     With a Lagrange multiplier Phi this becomes one linear elliptic equation (Sherman 1978):
//         d/dx(Rx dPhi/dx) + d/dy(Ry dPhi/dy) + d/dz(Rz dPhi/dz) = -(du0/dx + dv0/dy + dw0/dz)
//         Rx = Ry = 1/2,  Rz = alpha^2 / 2,  alpha = alpha_H / alpha_V (the stability ratio, 1 = neutral)
//     and the corrected wind is u = u0 + Rx dPhi/dx, v = v0 + Ry dPhi/dy, w = w0 + Rz dPhi/dz.
//  4. Solve. Finite elements on 8-node (trilinear) hexahedra with 2 x 2 x 2 Gauss quadrature. The mesh is
//     structured, so each node couples to at most its 27 neighbours and the matrix is stored as 27
//     coefficients per node. Boundary conditions as in WindNinja: Phi = 0 on the four sides and the top (air
//     flows in and out freely, the first guess is kept there); the ground is left natural (no flow through the
//     ground). The symmetric system is solved with conjugate gradients and a symmetric Gauss-Seidel (SSOR,
//     omega = 1) preconditioner, as WindNinja does.
//  5. Wind at the nodes. dPhi/dx etc. are discontinuous between elements, so each node takes the
//     inverse-distance-weighted average of the gradients at the surrounding Gauss points (WindNinja's
//     "stress smoothing", Thompson's book p. 228), then step 3's formulas give u, v, w.
//  6. 20-ft wind per fire cell. In each column the horizontal wind is read 6.1 m above the vegetation (scaled
//     down the log profile inside the lowest mesh cell; above it, interpolated along the log profile between
//     the two nodes either side), then interpolated bilinearly from the column centres to every fire cell.
//
// What this gives: speed-up over ridges and summits, slower wind in sheltered hollows, flow turned along
// valleys and round hills. What it cannot give (no momentum equation): flow separation, lee-side eddies and
// reversed flow; lee-side speeds come out too high. The published evaluations found the mass-conserving
// solver as good as a CFD solver at ridgetops and worse on lee slopes (Wagenbrenner et al. 2019,
// Atmosphere 10: 672). The wind is solved once, for one weather scenario.
// ============================================================================

struct WindSolverSettings {
    float columnSizeM = WIND_COLUMN_SIZE_M;   // Horizontal spacing of the wind mesh columns, m
    int layers = WIND_LAYERS;                 // Mesh nodes per column, ground to top (>= 3)
    float layerGrowth = 1.3f;                 // Each layer this much thicker than the one below (WindNinja: 1.3)
    float topAboveSummitM = WIND_TOP_ABOVE_SUMMIT_M;   // Domain top, at least this far above the highest column
    float maxAspectRatio = 400.0f;            // WindNinja's top: (column width / this) grown over the layers; used if higher
    float alpha = WIND_STABILITY_ALPHA;       // Stability ratio alpha_H / alpha_V: 1 neutral, > 1 unstable, < 1 stable
    float coriolisPerS = 1.0e-4f;             // |f| for the neutral boundary-layer height 0.2 u* / f (WindNinja's default when the latitude is unknown)
    double tolerance = 1.0e-4;                // CG stops when ||residual|| / ||right-hand side|| falls below this
    int maxIterations = 20000;                // Safety cap on CG iterations
};

// The solved wind on the fire grid, one entry per fire cell (row-major, GRID_WIDTH x GRID_HEIGHT)
struct WindField {
    std::vector<Vec2> toward;      // Unit vector the 20-ft wind blows TOWARDS (grid frame: x = east, y = south)
    std::vector<float> speedKmh;   // 20-ft wind speed (6.1 m above the vegetation), km/h

    // Diagnostics for the start-up log and the tests
    int columnsX = 0, columnsY = 0, layers = 0;   // Wind mesh size, nodes
    int iterations = 0;                           // CG iterations used
    double relResidual = 0.0;                     // Final ||residual|| / ||right-hand side|| (0 if nothing to correct)
    bool converged = false;                       // relResidual reached the tolerance
    double divergenceBefore = 0.0;                // RMS of div(u) of the first guess, 1/s (from the 3rd element layer up)
    double divergenceAfter = 0.0;                 // RMS of div(u) of the solved wind, 1/s (same elements)
};

// Solves the terrain wind for an open 20-ft wind of wind20Kmh blowing FROM windFromDeg (meteorological,
// clockwise from north). grid: GRID_WIDTH x GRID_HEIGHT cells with elevation, fuel model and canopy height
// filled in; cellSizeX/Y: ground size of one fire cell in metres.
WindField solveTerrainWind(const std::vector<Cell>& grid, float cellSizeX, float cellSizeY,
                           float wind20Kmh, float windFromDeg, const WindSolverSettings& settings = WindSolverSettings());
