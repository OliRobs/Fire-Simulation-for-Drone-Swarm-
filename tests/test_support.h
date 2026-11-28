#pragma once
// ============================================================================
// Helpers shared by the test files (tests/test_*.cpp) and, for the terrain generator, by tools/frame_bench.cpp:
// grid indexing, synthetic fires, a planned sweep for one fire, synthetic terrains and wind-test settings.
// Nothing in here uses the test framework, so the benchmark can link it too.
// ============================================================================
#include <cstdint>
#include <vector>
#include "config.h"
#include "distance_transform.h"
#include "fire_model.h"
#include "grid.h"
#include "mission_planner.h"
#include "perimeter.h"
#include "safe_planner.h"
#include "types.h"
#include "wind_solver.h"

// ----------------------------------------------------------------------------
// Grid indexing
// ----------------------------------------------------------------------------
// The grid most tests run on: the app's default size (config.h) with 10 m cells. Tests of grid handling also
// build their own small and non-square GridSpecs.
inline const GridSpec G{ GRID_WIDTH, GRID_HEIGHT };

// The same grid size with cells of cellX x cellY metres (tests that run the fire or the wind on 30 m cells)
inline GridSpec gridWithCells(float cellX, float cellY) { return GridSpec{ GRID_WIDTH, GRID_HEIGHT, cellX, cellY }; }
inline GridSpec gridWithCells(float cell) { return gridWithCells(cell, cell); }

// Cells on the grid G, and the index of cell (x, y)
constexpr int N = GRID_WIDTH * GRID_HEIGHT;
inline int at(int x, int y) { return y * GRID_WIDTH + x; }

// Centre of the grid. Tests with large fires (whose search area reaches ~2 km past the fire) put them here
// rather than at fixed cells, so they fit on any grid size in config.h, square or not (e.g. the current
// 1200 x 1000 cells).
constexpr int CX = GRID_WIDTH / 2;
constexpr int CY = GRID_HEIGHT / 2;

// ----------------------------------------------------------------------------
// Synthetic fires
// ----------------------------------------------------------------------------
// A filled disc of fire cells of radius r round (cx, cy), on grid G or on `spec`
std::vector<std::uint8_t> discFire(int cx, int cy, int r);
std::vector<std::uint8_t> discFire(const GridSpec& spec, int cx, int cy, int r);

// True if the segment a -> b passes through any fire cell (sampled every 0.1 cells), on grid G or on `spec`
bool crossesFire(const std::vector<std::uint8_t>& fire, Vec2 a, Vec2 b);
bool crossesFire(const GridSpec& spec, const std::vector<std::uint8_t>& fire, Vec2 a, Vec2 b);

// A fire with a concave bay: a disc with a bite taken out of its east side
std::vector<std::uint8_t> bayFire();

// Fires for the planner exactness tests: the bay (gap fill needs several passes), a fire with two bays, a fire
// on the grid's west edge (its perimeter line runs off the map and closes along the edge), a random blobby
// fire (many small bays and patches) and two separate fires (two loops)
std::vector<std::vector<std::uint8_t>> exactnessFires();

// ----------------------------------------------------------------------------
// Planning
// ----------------------------------------------------------------------------
// Everything a sweep test needs for one fire: the distance fields, the planner and the plan
struct SweepSetup {
    DistanceTransform live{ G };
    Perimeter perim;
    SafePlanner sp;
    MissionPlan plan;
    explicit SweepSetup(const std::vector<std::uint8_t>& fire);
};

// True if two lane lists are identical: same lanes in the same order, every point bit for bit
bool sameLanes(const std::vector<SweepLane>& a, const std::vector<SweepLane>& b);

// ----------------------------------------------------------------------------
// Fire spread
// ----------------------------------------------------------------------------
// The default weather with no wind
FireWeather calm();

// Furthest cell burned along a row/column from (x0, y0), in cells
int reach(const FireSimulation& sim, int x0, int y0, int dx, int dy);

// ----------------------------------------------------------------------------
// Synthetic terrain
// ----------------------------------------------------------------------------
// One generator for the ridge-and-valley maps the wind tests and the wind benchmark (wf_bench --wind) use.
// Elevation at ground position (X, Y) metres from the north-west corner:
//   z = base + ridgeAmp sin(X / ridgeLenX) cos(Y / ridgeLenY)      ridge-and-valley grain
//            + spurAmp sin((X + spurSlant Y) / spurLen)            smaller oblique spurs
//            + sum over peaks of height exp(-r^2 / (2 spread^2))   summits
// added up in that order. Optionally a forest patch (60 % cover, 18 m canopy) covers the north-west quarter,
// which gives the wind solver mixed roughness. Every other cell is the default (GR2 grass, no canopy).
struct Peak {
    double x, y;        // Centre, metres from the north-west corner
    double height;      // m
    double spread;      // Standard deviation, m
};
struct RidgeTerrain {
    double base = 0.0;
    double ridgeAmp = 0.0, ridgeLenX = 1.0, ridgeLenY = 1.0;
    double spurAmp = 0.0, spurSlant = 0.0, spurLen = 1.0;
    std::vector<Peak> peaks;
    bool forestNorthWest = false;
};
std::vector<Cell> ridgeValleyTerrain(const RidgeTerrain& t, float cellX, float cellY);

// The wind benchmark's map (about 1 km of relief and three summits; a San Gabriel-like scale)
RidgeTerrain windBenchTerrain();

// The wind tests' rugged map: about 700 m of relief, two summits placed at fixed fractions of the grid, and
// the forest patch, on cells of cellX x cellY metres
std::vector<Cell> ruggedWindTerrain(float cellX, float cellY);

// ----------------------------------------------------------------------------
// Terrain wind tests
// ----------------------------------------------------------------------------
// The wind tests run on 30 m cells with 600 m wind columns (60 x 50 x 20 nodes on the default grid) so each
// solve takes a fraction of a second. The same hill gives the same crest speed-up within ~2 % at the app's
// 10 m / 100 m.
constexpr float WIND_TEST_CELL_M = 30.0f;

// Solver settings for the wind tests: 600 m columns and the given stability alpha
WindSolverSettings windTestSettings(float alpha = 1.0f);

// A round Gaussian hill, 400 m high with a 2 km standard deviation, in the middle of the grid (30 m cells)
std::vector<Cell> gaussianHill();

// Compass heading (degrees clockwise from north) that a grid-frame vector points to
double headingDeg(Vec2 v);
