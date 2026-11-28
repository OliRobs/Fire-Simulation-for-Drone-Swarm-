// ============================================================================
// Helpers shared by the test files and the wind benchmark (see test_support.h)
// ============================================================================
#include "test_support.h"
#include <algorithm>
#include <cmath>
#include <random>

// ----------------------------------------------------------------------------
// Synthetic fires
// ----------------------------------------------------------------------------
std::vector<std::uint8_t> discFire(int cx, int cy, int r) {
    std::vector<std::uint8_t> m(N, 0);
    for (int y = 0; y < GRID_HEIGHT; ++y)
        for (int x = 0; x < GRID_WIDTH; ++x)
            if ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= r * r) m[at(x, y)] = 1;
    return m;
}

bool crossesFire(const std::vector<std::uint8_t>& fire, Vec2 a, Vec2 b) {
    const float len = std::sqrt(sqDist(a, b));
    const int steps = std::max(1, static_cast<int>(len / 0.1f));
    for (int i = 0; i <= steps; ++i) {
        const Vec2 p = a + (b - a) * (static_cast<float>(i) / steps);
        const int x = static_cast<int>(std::floor(p.x)), y = static_cast<int>(std::floor(p.y));
        if (x >= 0 && x < GRID_WIDTH && y >= 0 && y < GRID_HEIGHT && fire[at(x, y)]) return true;
    }
    return false;
}

std::vector<std::uint8_t> bayFire() {
    std::vector<std::uint8_t> fire = discFire(500, 500, 90);
    const auto bite = discFire(575, 500, 60);
    for (int i = 0; i < N; ++i) if (bite[i]) fire[i] = 0;
    return fire;
}

std::vector<std::vector<std::uint8_t>> exactnessFires() {
    std::vector<std::vector<std::uint8_t>> fires;
    fires.push_back(bayFire());
    std::vector<std::uint8_t> twoBays = discFire(CX, CY, 100);
    for (const auto& bite : { discFire(CX + 80, CY, 55), discFire(CX - 80, CY + 20, 50) })
        for (int i = 0; i < N; ++i) if (bite[i]) twoBays[i] = 0;
    fires.push_back(twoBays);
    fires.push_back(discFire(20, CY, 40));
    std::vector<std::uint8_t> blob(N, 0);
    std::mt19937 rng(77);
    std::uniform_int_distribution<int> off(-90, 90), rad(8, 40);
    for (int k = 0; k < 25; ++k) {
        const auto d = discFire(CX + off(rng), CY + off(rng), rad(rng));
        for (int i = 0; i < N; ++i) blob[i] |= d[i];
    }
    fires.push_back(blob);
    std::vector<std::uint8_t> twoFires = discFire(GRID_WIDTH / 4, GRID_HEIGHT / 4, 30);   // Two separate loops
    const auto other = discFire(GRID_WIDTH * 3 / 4, GRID_HEIGHT * 3 / 4, 35);
    for (int i = 0; i < N; ++i) twoFires[i] |= other[i];
    fires.push_back(twoFires);
    return fires;
}

// ----------------------------------------------------------------------------
// Planning
// ----------------------------------------------------------------------------
SweepSetup::SweepSetup(const std::vector<std::uint8_t>& fire) : perim(fire), sp(live) {
    live.compute(fire);
    std::mt19937 rng(2024);
    plan = planMission(live, perim.distance(), sp, rng);
}

bool sameLanes(const std::vector<SweepLane>& a, const std::vector<SweepLane>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i].pts.size() != b[i].pts.size() || a[i].loop != b[i].loop || a[i].gapFill != b[i].gapFill
            || a[i].seed.x != b[i].seed.x || a[i].seed.y != b[i].seed.y) return false;
        for (size_t k = 0; k < a[i].pts.size(); ++k)
            if (a[i].pts[k].x != b[i].pts[k].x || a[i].pts[k].y != b[i].pts[k].y) return false;
    }
    return true;
}

// ----------------------------------------------------------------------------
// Fire spread
// ----------------------------------------------------------------------------
FireWeather calm() {
    FireWeather w;
    w.wind20ftKmh = 0.0f;
    return w;
}

int reach(const FireSimulation& sim, int x0, int y0, int dx, int dy) {
    int k = 0;
    while (true) {
        const int x = x0 + (k + 1) * dx, y = y0 + (k + 1) * dy;
        if (x < 0 || x >= GRID_WIDTH || y < 0 || y >= GRID_HEIGHT || !sim.fireMask()[at(x, y)]) return k;
        ++k;
    }
}

// ----------------------------------------------------------------------------
// Synthetic terrain
// ----------------------------------------------------------------------------
std::vector<Cell> ridgeValleyTerrain(const RidgeTerrain& t, float cellX, float cellY) {
    std::vector<Cell> grid(N);
    for (int y = 0; y < GRID_HEIGHT; ++y)
        for (int x = 0; x < GRID_WIDTH; ++x) {
            const double X = x * cellX, Y = y * cellY;   // Metres from the north-west corner
            double z = t.base + t.ridgeAmp * std::sin(X / t.ridgeLenX) * std::cos(Y / t.ridgeLenY)
                     + t.spurAmp * std::sin((X + t.spurSlant * Y) / t.spurLen);
            for (const Peak& p : t.peaks) {
                const double dx = X - p.x, dy = Y - p.y;
                z += p.height * std::exp(-(dx * dx + dy * dy) / (2.0 * p.spread * p.spread));
            }
            Cell& c = grid[at(x, y)];
            c.elevation = static_cast<float>(z);
            if (t.forestNorthWest && x < GRID_WIDTH / 2 && y < GRID_HEIGHT / 2) {
                c.canopyCoverPct = 60.0f;
                c.canopyHeightM = 18.0f;
            }
        }
    return grid;
}

RidgeTerrain windBenchTerrain() {
    RidgeTerrain t;
    t.base = 900.0;
    t.ridgeAmp = 300.0;
    t.ridgeLenX = 1700.0;
    t.ridgeLenY = 1300.0;
    t.spurAmp = 150.0;
    t.spurSlant = 0.6;
    t.spurLen = 650.0;
    t.peaks = { { 3500, 2800, 450, 900 }, { 8200, 6100, 380, 1200 }, { 10400, 2200, 300, 700 } };
    return t;
}

std::vector<Cell> ruggedWindTerrain(float cellX, float cellY) {
    RidgeTerrain t;
    t.base = 500.0;
    t.ridgeAmp = 200.0;
    t.ridgeLenX = 2500.0;
    t.ridgeLenY = 1900.0;
    t.spurAmp = 80.0;
    t.spurSlant = 0.7;
    t.spurLen = 900.0;
    t.peaks = { { 0.3 * GRID_WIDTH * cellX, 0.4 * GRID_HEIGHT * cellY, 350.0, 2500.0 },
                { 0.7 * GRID_WIDTH * cellX, 0.6 * GRID_HEIGHT * cellY, 250.0, 1800.0 } };
    t.forestNorthWest = true;
    return ridgeValleyTerrain(t, cellX, cellY);
}

// ----------------------------------------------------------------------------
// Terrain wind tests
// ----------------------------------------------------------------------------
WindSolverSettings windTestSettings(float alpha) {
    WindSolverSettings s;
    s.columnSizeM = 600.0f;
    s.alpha = alpha;
    return s;
}

std::vector<Cell> gaussianHill() {
    std::vector<Cell> grid(N);
    const double sigmaCells = 2000.0 / WIND_TEST_CELL_M;
    for (int y = 0; y < GRID_HEIGHT; ++y)
        for (int x = 0; x < GRID_WIDTH; ++x) {
            const double r2 = ((x - CX) * (x - CX) + (y - CY) * (y - CY)) / (sigmaCells * sigmaCells);
            grid[at(x, y)].elevation = static_cast<float>(400.0 * std::exp(-0.5 * r2));
        }
    return grid;
}

double headingDeg(Vec2 v) { return std::atan2(v.x, -v.y) * 180.0 / 3.14159265358979; }
