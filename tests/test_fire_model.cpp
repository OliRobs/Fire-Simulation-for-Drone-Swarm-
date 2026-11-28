// ============================================================================
// Tests: Fire model (fire_model.h): spread on the grid, fuel lookup, wind convention, terrain from elevation and the 1 m DEM
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <initializer_list>
#include <string>
#include <vector>

#include "fire_model.h"
#include "fuel_models.h"
#include "rothermel.h"
#include "topography.h"
#include "viewport.h"
#include "test_framework.h"
#include "test_support.h"

// The fire simulation's bounding box must always be exactly the box of its fire mask
TEST(fire_bounds_track_the_fire_mask) {
    std::vector<Cell> grid(N);
    FireSimulation sim(G, grid);
    CHECK(sim.fireBounds().empty());
    sim.ignite(at(300, 400));
    for (int t = 0; t < 900; ++t) {
        sim.step();
        if (t == 300) sim.ignite(at(900, 120));   // A second fire, as a mouse click would start
        if (t % 150 != 0) continue;
        int x0, y0, x1, y1;
        CHECK(maskBounds(sim.fireMask(), GRID_WIDTH, GRID_HEIGHT, x0, y0, x1, y1));
        const GridRect b = sim.fireBounds();
        CHECK(b.x0 == x0 && b.y0 == y0 && b.x1 == x1 + 1 && b.y1 == y1 + 1);
    }
}

TEST(fuel_lookup) {
    for (int code : { 91, 92, 93, 98, 99, 0, -9999 }) {
        CHECK(isNonburnableCode(code));
        CHECK(findFuelModel(code) == nullptr);
    }
    CHECK(findFuelModel(150) == nullptr);           // Reserved, unused number: not guessed
    CHECK(std::string(findFuelModel(109)->name) == "GR9");
    CHECK(findFuelModel(106)->heatBtuLb == 9000);    // GR6 is the one model with 9000 Btu/lb
}

TEST(wind_direction_convention) {
    const Vec2 sw = windTowardGrid(225.0f);         // From the SW: blows to the NE (east = +x, north = -y)
    CHECK(std::fabs(sw.x - 0.7071f) < 1e-3f && std::fabs(sw.y + 0.7071f) < 1e-3f);
    const Vec2 w = windTowardGrid(270.0f);          // From the west: blows east
    CHECK(std::fabs(w.x - 1.0f) < 1e-4f && std::fabs(w.y) < 1e-4f);
}

TEST(topography_from_elevation) {
    // A plane rising 1 m per metre of northing (y decreasing): 45 degrees, uphill = -y
    std::vector<Cell> grid(N);
    for (int y = 0; y < GRID_HEIGHT; ++y)
        for (int x = 0; x < GRID_WIDTH; ++x) grid[at(x, y)].elevation = (GRID_HEIGHT - y) * 30.0f;
    computeTopography(gridWithCells(30.0f), grid, false);
    const Cell& c = grid[at(250, 250)];
    CHECK(std::fabs(c.slopeDeg - 45.0f) < 0.01f);
    CHECK(std::fabs(c.upslopeX) < 1e-5f && std::fabs(c.upslopeY + 1.0f) < 1e-5f);
}

// --- Simulation-cell terrain from a finer (1 m) DEM ---
namespace {
    // Fine DEM buffer for cellsX x cellsY cells of `factor` pixels, with the one-pixel halo, filled from
    // z(px, py) where px, py are pixel columns/rows counted from the halo's top-left corner
    std::vector<float> fineDem(int cellsX, int cellsY, int factor, const std::function<float(int, int)>& z) {
        const int w = cellsX * factor + 2, h = cellsY * factor + 2;
        std::vector<float> buf(static_cast<size_t>(w) * h);
        for (int py = 0; py < h; ++py)
            for (int px = 0; px < w; ++px) buf[static_cast<size_t>(py) * w + px] = z(px, py);
        return buf;
    }
}

// A tilted plane: every 1 m pixel has the same slope, so each 10 m cell must report exactly that slope, the
// plane's uphill direction and the mean elevation of its 100 pixels
TEST(fine_dem_plane_gives_exact_cell_terrain) {
    // Rises 0.1 m per metre east and 0.2 m per metre north (rows run south, so -0.2 per row)
    auto plane = [](int px, int py) { return 500.0f + 0.1f * px - 0.2f * py; };
    std::vector<CellTopography> out;
    aggregateFineTopography(fineDem(3, 2, 10, plane), 3, 2, 10, 1.0f, 1.0f, out);
    CHECK(out.size() == 6);
    const float expectSlope = static_cast<float>(std::atan(std::sqrt(0.01 + 0.04)) * 57.29577951308232);  // 12.6 deg
    for (int k = 0; k < 6; ++k) {
        CHECK(out[k].validPixels == 100);
        CHECK(std::fabs(out[k].slopeDeg - expectSlope) < 0.01f);
        CHECK(std::fabs(out[k].upslopeX - 0.4472136f) < 1e-4f && std::fabs(out[k].upslopeY + 0.8944272f) < 1e-4f);
    }
    // Cell (1, 1): pixels 11..20 across (in buffer columns, after the halo) and 11..20 down, centre 15.5
    CHECK(std::fabs(out[4].elevation - plane(0, 0) - 0.1f * 15.5f + 0.2f * 15.5f) < 1e-3f);
}

// Corrugated ground (1 m ridges and gullies at 45 degrees) averages to a FLAT 10 m surface, but the fire
// still sees steep slopes: the cell slope must come from the 1 m pixels, not from the smoothed elevations
TEST(fine_dem_keeps_steep_detail) {
    auto ridges = [](int px, int) { return 100.0f + std::fabs(static_cast<float>(px % 4) - 2.0f); };
    std::vector<CellTopography> out;
    aggregateFineTopography(fineDem(2, 2, 10, ridges), 2, 2, 10, 1.0f, 1.0f, out);
    // Horn's 3x3 difference gives 45 deg on the sides and 0 on each ridge top and gully floor: 22.5 deg mean
    for (const CellTopography& t : out) CHECK(std::fabs(t.slopeDeg - 22.5f) < 0.01f);

    // The same ground at 10 m: cell means are nearly equal, so a 10 m slope would be close to flat
    std::vector<Cell> coarse(N);
    for (int i = 0; i < N; ++i) coarse[i].elevation = 101.0f;
    for (int k = 0; k < 4; ++k) coarse[at(k % 2, k / 2)].elevation = out[k].elevation;
    computeTopography(gridWithCells(10.0f), coarse, false);
    CHECK(coarse[at(0, 0)].slopeDeg < 5.0f);
}

// Missing DEM pixels (NoData) are left out of the averages, and a cell with none at all is marked uncovered
TEST(fine_dem_skips_missing_pixels) {
    auto flat = [](int px, int) { return px < 11 ? NAN : 250.0f; };   // Left cell (buffer columns 0-10) has no data
    std::vector<CellTopography> out;
    aggregateFineTopography(fineDem(2, 1, 10, flat), 2, 1, 10, 1.0f, 1.0f, out);
    CHECK(out[0].validPixels == 0);
    CHECK(out[1].validPixels == 100);
    CHECK(std::fabs(out[1].elevation - 250.0f) < 1e-4f);
    CHECK(out[1].slopeDeg < 1e-3f);   // A NaN neighbour is replaced by the centre, so no false cliff at the edge
}

// --- Fire spread on the grid ---

// No wind, flat, uniform GR9: the front must advance at R_0 (along rows and along diagonals)
TEST(fire_front_moves_at_rothermel_rate) {
    std::vector<Cell> grid(N);
    for (Cell& c : grid) c.fuelModel = 109;
    FireSimulation sim(gridWithCells(30.0f), grid, calm());
    CHECK(sim.ignite(at(250, 250)));
    const int ticks = 1200;                          // 200 minutes
    for (int t = 0; t < ticks; ++t) sim.step();

    const double r0 = computeFuelBed(*findFuelModel(109), FuelMoisture()).r0 * FT_PER_MIN_TO_M_PER_MIN;
    const double expectedCells = r0 * sim.minutes() / 30.0;
    const int east = reach(sim, 250, 250, 1, 0), diag = reach(sim, 250, 250, 1, 1);
    std::printf("  (expected %.1f cells, reached %d along a row, %.1f along a diagonal)\n",
                expectedCells, east, diag * std::sqrt(2.0));
    CHECK(std::abs(east - expectedCells) <= 1.0);
    CHECK(std::fabs(diag * std::sqrt(2.0) - expectedCells) <= 1.5);
}

TEST(fire_spreads_downwind_and_upslope) {
    // Wind from the west: the fire runs east
    {
        std::vector<Cell> grid(N);
        FireWeather w;
        w.windFromDeg = 270.0f;
        FireSimulation sim(gridWithCells(30.0f), grid, w);
        sim.ignite(at(250, 250));
        for (int t = 0; t < 1500; ++t) sim.step();
        CHECK(reach(sim, 250, 250, 1, 0) > 3 * std::max(1, reach(sim, 250, 250, -1, 0)));
    }
    // No wind, ground rising to the north: the fire runs uphill
    {
        std::vector<Cell> grid(N);
        for (int y = 0; y < GRID_HEIGHT; ++y)
            for (int x = 0; x < GRID_WIDTH; ++x) grid[at(x, y)].elevation = (GRID_HEIGHT - y) * 15.0f;
        computeTopography(gridWithCells(30.0f), grid, false);
        FireSimulation sim(gridWithCells(30.0f), grid, calm());
        sim.ignite(at(250, 250));
        for (int t = 0; t < 3000; ++t) sim.step();
        CHECK(reach(sim, 250, 250, 0, -1) > 2 * std::max(1, reach(sim, 250, 250, 0, 1)));
    }
}

TEST(fire_never_enters_nonburnable_cells) {
    std::vector<Cell> grid(N);
    for (Cell& c : grid) c.fuelModel = 104;          // GR4: fast grass
    // A vertical strip of water (NB8) down the middle of the map
    for (int y = 0; y < GRID_HEIGHT; ++y)
        for (int x = 240; x < 245; ++x) grid[at(x, y)].fuelModel = 98;

    FireWeather w;
    w.windFromDeg = 270.0f;                          // Wind blowing towards the water
    FireSimulation sim(G, grid, w);
    CHECK(sim.ignite(at(230, 250)));
    CHECK(!sim.ignite(at(242, 250)));               // Can't light water by clicking on it
    for (int t = 0; t < 3000; ++t) sim.step();

    CHECK(sim.fireMask()[at(239, 250)] == 1);       // It did reach the water...
    for (int y = 0; y < GRID_HEIGHT; ++y)
        for (int x = 240; x < GRID_WIDTH; ++x) CHECK(sim.grid()[at(x, y)].state == CellState::Unburned);
}

// A one-cell-wide DIAGONAL road must hold too: the fire may not squeeze between its corners
TEST(roads_are_fuel_breaks_even_diagonally) {
    std::vector<Cell> grid(N);
    for (Cell& c : grid) c.fuelModel = 104;
    // The road runs corner to corner across the largest square in the top-left of the map, from (0, S - 1) to
    // (S - 1, 0) with S the shorter side, so it meets two map edges and cuts the map in two whatever its shape
    const int S = std::min(GRID_WIDTH, GRID_HEIGHT);
    for (int i = 0; i < S; ++i) grid[at(i, S - 1 - i)].roadClass = 22;   // Tertiary road, SW-NE

    FireWeather w;
    w.windFromDeg = 315.0f;                          // From the NW, straight at the road
    FireSimulation sim(G, grid, w);
    // Positions relative to the map size, so the test holds for any grid in config.h:
    // k + k = S - 2 is the diagonal just north-west of the road
    const int k = (S - 2) / 2;
    CHECK(!sim.canBurn(at(k + 1, k)));               // On the road (x + y = S - 1)
    CHECK(sim.ignite(at(k - 49, k - 49)));           // North-west of the road, ~50 cells from it
    for (int t = 0; t < 4000; ++t) sim.step();
    int beyond = 0;
    for (int y = 0; y < GRID_HEIGHT; ++y)
        for (int x = 0; x < GRID_WIDTH; ++x)
            if (x + y > S - 1 && sim.fireMask()[at(x, y)]) ++beyond;
    CHECK(sim.fireMask()[at(k, k)] == 1);           // The fire reached the road...
    CHECK(beyond == 0);                              // ...and didn't cross it
}

TEST(fire_is_deterministic) {
    std::vector<Cell> grid(N);
    FireSimulation a(G, grid), b(G, grid);
    a.ignite(at(250, 250)); b.ignite(at(250, 250));
    for (int t = 0; t < 600; ++t) { a.step(); b.step(); }
    CHECK(a.fireMask() == b.fireMask());
    // The tick length only changes how often we look, not where the fire goes
    FireSimulation c(G, grid, FireWeather(), SIM_MINUTES_PER_TICK * 4);
    c.ignite(at(250, 250));
    for (int t = 0; t < 150; ++t) c.step();
    CHECK(c.fireMask() == a.fireMask());
}
