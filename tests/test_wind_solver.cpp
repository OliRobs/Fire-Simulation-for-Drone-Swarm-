// ============================================================================
// Tests: Terrain-adjusted wind (wind_solver.h)
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

#include "fire_model.h"
#include "wind_solver.h"
#include "test_framework.h"
#include "test_support.h"

// ============================================================================
// Terrain-adjusted wind (wind_solver.h)
// ============================================================================
// These run on 30 m cells with 600 m wind columns (60 x 50 x 20 nodes on the default grid) so each solve takes
// a fraction of a second. The same hill gives the same crest speed-up within ~2 % at the app's 10 m / 100 m.

// Flat ground with even fuel: the open wind is already mass-conserving, so nothing may change
TEST(terrain_wind_keeps_the_open_wind_on_flat_ground) {
    std::vector<Cell> grid(N);   // Flat, GR2 everywhere
    const WindField f = solveTerrainWind(gridWithCells(WIND_TEST_CELL_M), grid, 20.0f, 225.0f, windTestSettings());
    const Vec2 open = windTowardGrid(225.0f);
    double worstSpeed = 0.0, worstDir = 0.0;
    for (int i = 0; i < N; ++i) {
        worstSpeed = std::max(worstSpeed, std::fabs(f.speedKmh[i] - 20.0));
        worstDir = std::max(worstDir, static_cast<double>(std::fabs(f.toward[i].x - open.x) + std::fabs(f.toward[i].y - open.y)));
    }
    std::printf("  (%d iterations, worst speed error %.2e km/h, worst direction error %.2e)\n", f.iterations, worstSpeed, worstDir);
    CHECK(f.converged);
    CHECK(f.iterations == 0);   // Nothing to correct
    CHECK(worstSpeed < 1e-3);
    CHECK(worstDir < 1e-5);
}

// Wind from the west over a round hill: faster over the crest, slower at the upwind foot and in the lee, turned
// out round both flanks, mirror-symmetric north/south, and mass conserved far better than the first guess
TEST(terrain_wind_speeds_up_over_a_hill_and_turns_round_it) {
    const WindField f = solveTerrainWind(gridWithCells(WIND_TEST_CELL_M), gaussianHill(), 20.0f, 270.0f, windTestSettings());
    const int d = static_cast<int>(2.0 * 2000.0 / WIND_TEST_CELL_M);   // Two standard deviations, in cells
    const float crest = f.speedKmh[at(CX, CY)], upwind = f.speedKmh[at(CX - d, CY)], lee = f.speedKmh[at(CX + d, CY)];
    const Vec2 nw = f.toward[at(CX - d, CY - d)], sw = f.toward[at(CX - d, CY + d)];
    std::printf("  (%d iterations; crest %.2f, upwind foot %.2f, lee foot %.2f km/h; NW flank heading %.1f, SW %.1f deg;"
                " RMS div %.2e -> %.2e)\n", f.iterations, crest, upwind, lee, headingDeg(nw), headingDeg(sw),
                f.divergenceBefore, f.divergenceAfter);
    CHECK(f.converged);
    CHECK(crest > 1.15f * 20.0f);             // Speed-up over the crest (~19 %)
    CHECK(upwind < 0.95f * 20.0f);            // Slowed where the air meets the slope
    CHECK(lee < 0.95f * 20.0f);               // Mass-conserving: no wake, the lee mirrors the windward side
    CHECK(nw.y < -0.02f);                     // North flank: turned north (-y) round the hill
    CHECK(sw.y > 0.02f);                      // South flank: turned south
    // ...by the same amount, and the two flanks equally fast: the hill is symmetric. (Not to the last digit:
    // the hill is centred on a fire cell, half a cell off the middle of the 20-cell wind columns.)
    CHECK(std::fabs(nw.y + sw.y) < 1e-3f);
    CHECK(std::fabs(f.speedKmh[at(CX, CY - d)] - f.speedKmh[at(CX, CY + d)]) < 0.01f);
    CHECK(f.divergenceAfter < 0.1 * f.divergenceBefore);   // ~19 times smaller
}

// A north-south valley with a south-west wind: on the valley floor the wind slows and turns to run along the
// valley (towards heading 0), while the open wind heads 45 deg
TEST(terrain_wind_is_channelled_along_a_valley) {
    std::vector<Cell> grid(N);
    for (int y = 0; y < GRID_HEIGHT; ++y)
        for (int x = 0; x < GRID_WIDTH; ++x) {
            const double dx = (x - CX) * WIND_TEST_CELL_M;   // Metres from the valley axis
            grid[at(x, y)].elevation = static_cast<float>(600.0 * (1.0 - std::exp(-dx * dx / (2.0 * 1500.0 * 1500.0))));
        }
    const WindField f = solveTerrainWind(gridWithCells(WIND_TEST_CELL_M), grid, 20.0f, 225.0f, windTestSettings());
    const double floorHeading = headingDeg(f.toward[at(CX, CY)]);
    std::printf("  (valley floor: %.2f km/h heading %.1f deg)\n", f.speedKmh[at(CX, CY)], floorHeading);
    CHECK(f.converged);
    CHECK(floorHeading > 0.0 && floorHeading < 35.0);   // Turned at least 10 deg towards the valley axis
    CHECK(f.speedKmh[at(CX, CY)] < 0.9f * 20.0f);        // Sheltered floor
}

// Stability: stable air (alpha < 1) is pushed round the hill more than neutral air, unstable air (alpha > 1)
// goes over it and is turned less (WindNinja's convention, alpha = alpha_H / alpha_V)
SLOW_TEST(terrain_wind_stability_changes_how_much_air_goes_round) {
    const std::vector<Cell> hill = gaussianHill();
    const int d = static_cast<int>(2.0 * 2000.0 / WIND_TEST_CELL_M);
    double turn[3];
    const float alphas[3] = { 0.3f, 1.0f, 3.0f };
    for (int i = 0; i < 3; ++i) {
        const WindField f = solveTerrainWind(gridWithCells(WIND_TEST_CELL_M), hill, 20.0f, 270.0f, windTestSettings(alphas[i]));
        CHECK(f.converged);
        turn[i] = -f.toward[at(CX - d, CY - d)].y;   // How far north the NW flank wind is turned
    }
    std::printf("  (NW flank turn: stable %.3f, neutral %.3f, unstable %.3f)\n", turn[0], turn[1], turn[2]);
    CHECK(turn[0] > turn[1] && turn[1] > turn[2]);
}

// Golden results: six solves recorded before the reference kernels (the solver as first written) were removed,
// when the fast paths were proven to give the reference's results bit for bit on exactly these meshes. Each is
// a hash of the whole solved field (windFieldHash: iterations, residual, every cell's speed and direction), so
// any change to the solved wind, down to the last bit, fails here. The cases cover: square and non-square
// cells; meshes whose free rows are and are not a multiple of the sweep band (4 rows); meshes narrower than the
// band's skew (2 free nodes per row); 3 to 20 layers; stable, neutral and unstable air (alpha 0.3 / 1 / 3, so
// Rz is not always a power of two); and mixed roughness (a forest patch).
//
// The hashes depend on the compiler's floating point (and libm, which builds the terrain), so they are only
// compared where they were recorded: GCC on x86-64 Linux without FMA contraction, as in the Linux CI job (the
// same condition as tools/wind_checksum.txt). Elsewhere (MSVC CI) the solves still run and must converge.
// A change MEANT to alter the wind updates these values and tools/wind_checksum.txt together, and says why.
#if defined(__GNUC__) && !defined(__clang__) && defined(__x86_64__) && defined(__linux__) && !defined(__FMA__)
#define WF_COMPARE_WIND_GOLDEN 1
#else
#define WF_COMPARE_WIND_GOLDEN 0
#endif

SLOW_TEST(terrain_wind_matches_the_recorded_golden_fields) {
    struct Case { float cellX, cellY, columnM; int layers; float alpha; float fromDeg; bool hill; std::uint64_t golden; };
    const Case cases[] = {
        // 60 x 50 x 20 nodes: 48 free rows (a multiple of 4)
        { 30.0f, 30.0f, 600.0f, 20, 1.0f, 270.0f, true, 0x7602f37c60349b30ULL },
        // Stable air, rugged terrain with a forest patch
        { 30.0f, 30.0f, 600.0f, 20, 0.3f, 225.0f, false, 0x4de32e12c1bf6794ULL },
        // Non-square cells, unstable air: 48 x 34 x 7 (32 free rows)
        { 30.0f, 25.0f, 750.0f, 7, 3.0f, 200.0f, false, 0x224151796e88663eULL },
        // 24 x 20 x 5: 18 free rows (4 bands + 2)
        { 30.0f, 30.0f, 1500.0f, 5, 1.0f, 135.0f, false, 0xd1a40a20e3ee939cULL },
        // 6 x 5 x 4: 3 free rows (one partial band)
        { 30.0f, 30.0f, 6000.0f, 4, 1.0f, 45.0f, false, 0xbdd5fa1e54f6f02dULL },
        // 4 x 4 x 3: 2 free nodes per row, narrower than the skew
        { 30.0f, 30.0f, 9000.0f, 3, 1.0f, 315.0f, false, 0x26ba2bdd8b011e7bULL },
    };
    int checked = 0;
    for (const Case& c : cases) {
        const std::vector<Cell> grid = c.hill ? gaussianHill() : ruggedWindTerrain(c.cellX, c.cellY);
        WindSolverSettings s = windTestSettings(c.alpha);
        s.columnSizeM = c.columnM;
        s.layers = c.layers;
        const WindField f = solveTerrainWind(gridWithCells(c.cellX, c.cellY), grid, 20.0f, c.fromDeg, s);
        const std::uint64_t hash = windFieldHash(f);
        std::printf("  (%d x %d x %d mesh, alpha %.1f: %d CG iterations, hash %016llx%s)\n", f.columnsX, f.columnsY,
                    f.layers, c.alpha, f.iterations, static_cast<unsigned long long>(hash),
                    !WF_COMPARE_WIND_GOLDEN ? ", not compared on this compiler" : hash == c.golden ? "" : " DIFFERENT");
        CHECK(f.iterations > 0);   // The solve really ran (not the flat-ground shortcut)
        CHECK(f.converged);
        if (WF_COMPARE_WIND_GOLDEN) CHECK(hash == c.golden);
        ++checked;
    }
    CHECK(checked == 6);
}

// The fire model reads each cell's wind from the field: a field equal to the open wind changes nothing, and a
// field blowing north drives the fire north even though the weather's own wind blows east
TEST(fire_uses_the_terrain_wind_field) {
    std::vector<Cell> grid(N);
    FireWeather w;
    w.windFromDeg = 270.0f;   // The weather's open wind blows east

    // A field that is exactly the open wind: identical spread rates to running without a field
    WindField same;
    same.toward.assign(N, windTowardGrid(w.windFromDeg));
    same.speedKmh.assign(N, w.wind20ftKmh);
    const FireSimulation plain(gridWithCells(WIND_TEST_CELL_M), grid, w);
    const FireSimulation withSame(gridWithCells(WIND_TEST_CELL_M), grid, w, SIM_MINUTES_PER_TICK, &same);
    bool identical = true;
    for (int i = 0; i < N && identical; i += 997)
        for (int dir = 0; dir < NUM_DIRS; ++dir) identical = identical && plain.directionalRos(i, dir) == withSame.directionalRos(i, dir);
    CHECK(identical);

    // A field blowing north: the fire runs north, not east
    WindField north;
    north.toward.assign(N, windTowardGrid(180.0f));   // From the south: blows north (-y)
    north.speedKmh.assign(N, w.wind20ftKmh);
    FireSimulation sim(gridWithCells(WIND_TEST_CELL_M), grid, w, SIM_MINUTES_PER_TICK, &north);
    sim.ignite(at(250, 250));
    for (int t = 0; t < 1500; ++t) sim.step();
    const int up = reach(sim, 250, 250, 0, -1), east = reach(sim, 250, 250, 1, 0);
    std::printf("  (north %d cells, east %d cells)\n", up, east);
    CHECK(up > 3 * std::max(1, east));
}
