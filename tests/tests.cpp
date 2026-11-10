// ============================================================================
// Headless tests for the simulation core (no SFML, no GDAL).
// Build: part of the CMake project (target wf_tests); run with ctest or ./wf_tests.
// Each TEST is a plain function; CHECK records a failure and carries on, so one run reports everything.
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <initializer_list>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "aco.h"
#include "distance_transform.h"
#include "drone.h"
#include "fire_model.h"
#include "crown_fire.h"
#include "fuel_models.h"
#include "rothermel.h"
#include "topography.h"
#include "viewport.h"
#include "mission_controller.h"
#include "mission_planner.h"
#include "perimeter.h"
#include "safe_planner.h"
#include "voronoi.h"

// --- Minimal test framework ---
namespace {
    int g_failures = 0;
    std::vector<std::pair<const char*, std::function<void()>>>& registry() {
        static std::vector<std::pair<const char*, std::function<void()>>> r;
        return r;
    }
    struct Registrar {
        Registrar(const char* name, std::function<void()> fn) { registry().emplace_back(name, std::move(fn)); }
    };
}
#define TEST(name) static void name(); static Registrar reg_##name(#name, name); static void name()
#define CHECK(cond) do { if (!(cond)) { ++g_failures; \
    std::printf("  FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while (0)

// --- Helpers ---
namespace {
    constexpr int N = GRID_WIDTH * GRID_HEIGHT;
    int at(int x, int y) { return y * GRID_WIDTH + x; }

    // A filled disc of fire cells
    std::vector<std::uint8_t> discFire(int cx, int cy, int r) {
        std::vector<std::uint8_t> m(N, 0);
        for (int y = 0; y < GRID_HEIGHT; ++y)
            for (int x = 0; x < GRID_WIDTH; ++x)
                if ((x - cx) * (x - cx) + (y - cy) * (y - cy) <= r * r) m[at(x, y)] = 1;
        return m;
    }

    // True if the segment a -> b passes through any fire cell (sampled every 0.1 cells)
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
}

// ============================================================================
// Distance transform
// ============================================================================
TEST(edt_matches_brute_force) {
    const int w = 37, h = 23;
    std::mt19937 rng(7);
    std::vector<std::uint8_t> src(w * h, 0);
    for (int k = 0; k < 12; ++k) src[rng() % (w * h)] = 1;

    DistanceTransform edt(w, h);
    edt.compute(src);
    for (int i = 0; i < w * h; ++i) {
        double best = 1e18;
        for (int j = 0; j < w * h; ++j) {
            if (!src[j]) continue;
            const double dx = i % w - j % w, dy = i / w - j / w;
            best = std::min(best, dx * dx + dy * dy);
        }
        CHECK(std::fabs(edt.squaredDistance(i) - best) < 1e-9);
    }
}

// ============================================================================
// Safe planner
// ============================================================================
TEST(planner_straight_line_when_clear) {
    DistanceTransform edt(GRID_WIDTH, GRID_HEIGHT);
    edt.compute(discFire(400, 400, 5));
    SafePlanner sp(edt);
    const auto p = sp.leg(Vec2(50.5f, 50.5f), Vec2(150.5f, 60.5f));
    CHECK(p.size() == 1);
    CHECK(!p.empty() && sqDist(p.back(), Vec2(150.5f, 60.5f)) < 1e-6f);
}

TEST(planner_detours_around_fire) {
    const auto fire = discFire(250, 250, 20);
    DistanceTransform edt(GRID_WIDTH, GRID_HEIGHT);
    edt.compute(fire);
    SafePlanner sp(edt);

    const Vec2 a(200.5f, 250.5f), b(300.5f, 250.5f);   // Straight line goes through the middle of the fire
    CHECK(!sp.lineClear(a, b));
    const auto p = sp.leg(a, b);
    CHECK(!p.empty());
    CHECK(!p.empty() && sqDist(p.back(), b) < 1e-6f);

    // Every segment of the detour must stay out of the hard no-fly zone
    Vec2 prev = a;
    for (const Vec2& q : p) { CHECK(sp.lineClear(prev, q)); prev = q; }
}

TEST(planner_reports_no_path_instead_of_straight_line) {
    // Goal enclosed by a ring of fire: the old code returned a straight line straight through the fire
    std::vector<std::uint8_t> ring(N, 0);
    for (int y = 0; y < GRID_HEIGHT; ++y)
        for (int x = 0; x < GRID_WIDTH; ++x) {
            const int r2 = (x - 250) * (x - 250) + (y - 250) * (y - 250);
            if (r2 >= 28 * 28 && r2 <= 32 * 32) ring[at(x, y)] = 1;
        }
    DistanceTransform edt(GRID_WIDTH, GRID_HEIGHT);
    edt.compute(ring);
    SafePlanner sp(edt);
    CHECK(sp.leg(Vec2(50.5f, 250.5f), Vec2(250.5f, 250.5f)).empty());

    // A second search after a failed one must still work (A* buffers are reused between calls)
    CHECK(!sp.leg(Vec2(50.5f, 250.5f), Vec2(450.5f, 250.5f)).empty());
}

TEST(planner_goal_inside_no_fly_zone_stops_short) {
    DistanceTransform edt(GRID_WIDTH, GRID_HEIGHT);
    edt.compute(discFire(250, 250, 10));
    SafePlanner sp(edt);
    const Vec2 goal(251.5f, 250.5f);                     // On the fire itself
    const auto p = sp.leg(Vec2(100.5f, 250.5f), goal);
    CHECK(!p.empty());
    if (!p.empty()) {
        const Vec2 end = p.back();
        CHECK(sp.cellOk(static_cast<int>(end.x), static_cast<int>(end.y)));  // Ends in flyable space
        CHECK(sqDist(end, goal) > 1.0f);
    }
}

// ============================================================================
// Fire model
// ============================================================================
// --- Rothermel model ---
namespace {
    // Scott & Burgan's adjective spread-rate classes (RMRS-GTR-153 Table 5), in chains per hour
    enum SpreadClass { VL, L, M, H, VH, E };
    SpreadClass classify(double chainsPerHour) {
        if (chainsPerHour < 2) return VL;
        if (chainsPerHour < 5) return L;
        if (chainsPerHour < 20) return M;
        if (chainsPerHour < 50) return H;
        if (chainsPerHour < 150) return VH;
        return E;
    }
    // Head fire rate (ft/min) for a fuel model on flat ground with a midflame wind blowing along +x
    double headRosFtMin(const FuelModel& m, const FuelMoisture& mo, double midflameMph, bool windLimit) {
        const SurfaceFuelBed bed = computeFuelBed(m, mo);
        return computeSpread(bed, midflameMph * FT_PER_MIN_PER_MPH, 1, 0, 0, 0, 0, windLimit).headRos;
    }
}

// Every fuel model's head fire rate must fall in the spread class the selection guide gives it
// (RMRS-GTR-153 pp. 10-12: moisture scenario D2L2, 5 mi/h midflame wind, zero slope).
TEST(rothermel_reproduces_scott_burgan_spread_classes) {
    const std::pair<const char*, SpreadClass> expected[] = {
        { "GR1", M }, { "GR2", H }, { "GR3", H }, { "GR4", VH }, { "GR5", VH }, { "GR6", VH }, { "GR7", VH },
        { "GR8", VH }, { "GR9", E }, { "GS1", M }, { "GS2", H }, { "GS3", H }, { "GS4", H },
        { "SH1", VL }, { "SH2", L }, { "SH3", L }, { "SH4", H }, { "SH5", VH }, { "SH6", H }, { "SH7", H },
        { "SH8", H }, { "SH9", H }, { "TU1", L }, { "TU2", M }, { "TU3", H }, { "TU4", M }, { "TU5", M },
        { "TL1", VL }, { "TL2", VL }, { "TL3", VL }, { "TL4", L }, { "TL5", L }, { "TL6", M }, { "TL7", L },
        { "TL8", M }, { "TL9", M }, { "SB1", M }, { "SB2", M }, { "SB3", H }, { "SB4", VH },
    };
    int count = 0;
    const FuelModel* models = allFuelModels(count);
    CHECK(count == 40);
    for (int i = 0; i < count; ++i) {
        for (const auto& [name, cls] : expected) {
            if (std::string(name) != models[i].name) continue;
            const double chainsPerHour = headRosFtMin(models[i], FuelMoisture(), 5.0, true) * 60.0 / 66.0;
            if (classify(chainsPerHour) != cls) {
                std::printf("  %s: %.1f ch/h\n", name, chainsPerHour);
                CHECK(classify(chainsPerHour) == cls);
            }
        }
    }
}

// Packing ratio and characteristic SAV against the summary values on the fuel model pages (RMRS-GTR-153)
TEST(rothermel_fuel_bed_matches_published_values) {
    struct Ref { int code; double beta; double sigma; };
    for (const Ref r : { Ref{ 101, 0.00143, 2054 }, Ref{ 142, 0.01198, 1672 }, Ref{ 181, 0.04878, 1716 },
                         Ref{ 164, 0.01865, 2216 }, Ref{ 146, 0.00412, 1144 } }) {
        const SurfaceFuelBed bed = computeFuelBed(*findFuelModel(r.code), FuelMoisture());
        CHECK(std::fabs(bed.packingRatio - r.beta) / r.beta < 0.005);
        CHECK(std::fabs(bed.sigma - r.sigma) < 1.0);
    }
}

TEST(rothermel_moisture_and_curing) {
    const FuelModel& gr2 = *findFuelModel(102);
    // Dead fuel wetter than the 15 % moisture of extinction: no spread at all
    FuelMoisture wet;
    wet.dead1h = wet.dead10h = wet.dead100h = 0.30;
    CHECK(!computeFuelBed(gr2, wet).burnable());
    // Dynamic model: fully cured grass (herb moisture 30 %) burns much faster than green grass (120 %)
    FuelMoisture cured, green;
    cured.liveHerb = 0.30;
    green.liveHerb = 1.20;
    CHECK(computeFuelBed(gr2, cured).r0 > 2.0 * computeFuelBed(gr2, green).r0);
}

TEST(rothermel_wind_slope_and_ellipse) {
    const SurfaceFuelBed bed = computeFuelBed(*findFuelModel(145), FuelMoisture());   // SH5
    // No wind, no slope: a circle, spreading at R_0 every way
    const SpreadEllipse still = computeSpread(bed, 0, 1, 0, 0, 0, 0, false);
    CHECK(std::fabs(still.headRos - bed.r0) < 1e-9 && still.eccentricity == 0.0);
    CHECK(std::fabs(still.rosToward(0, 1) - bed.r0) < 1e-9);

    // Wind towards +x: head fire faster, flanks slower, backing slowest
    const SpreadEllipse windy = computeSpread(bed, 5 * FT_PER_MIN_PER_MPH, 1, 0, 0, 0, 0, false);
    CHECK(windy.headRos > 3.0 * bed.r0);
    CHECK(std::fabs(windy.rosToward(1, 0) - windy.headRos) < 1e-9);
    CHECK(windy.rosToward(0, 1) < windy.headRos && windy.backingRos() < windy.rosToward(0, 1));

    // Slope alone (30 degrees, uphill = -y): head points uphill and gets the slope factor
    const double tan30 = std::tan(30.0 * 3.14159265358979 / 180.0);
    const SpreadEllipse slope = computeSpread(bed, 0, 1, 0, tan30, 0, -1, false);
    CHECK(std::fabs(slope.dirY + 1.0) < 1e-9);
    CHECK(std::fabs(slope.headRos - bed.r0 * (1 + bed.slopeFactor(tan30))) < 1e-9);

    // Wind across the slope: the head direction lies between the two
    const SpreadEllipse both = computeSpread(bed, 5 * FT_PER_MIN_PER_MPH, 1, 0, tan30, 0, -1, false);
    CHECK(both.dirX > 0.0 && both.dirY < 0.0);

    // Unsheltered WAF: deeper fuel beds feel more of the 20-ft wind (Andrews 2012)
    CHECK(unshelteredWaf(1.0) > 0.3 && unshelteredWaf(1.0) < 0.4);
    CHECK(unshelteredWaf(6.0) > unshelteredWaf(1.0));
}

// --- Canopy and crown fire ---
TEST(fuel_model_10_matches_published_values) {
    // RMRS-GTR-371 Table 8a: characteristic SAV 1764 1/ft, relative packing ratio 2.35
    const SurfaceFuelBed fm10 = computeFuelBed(fuelModel10(), FuelMoisture());
    CHECK(std::fabs(fm10.sigma - 1764.0) < 1.0);
    CHECK(std::fabs(fm10.relativePacking - 2.35) < 0.01);
}

TEST(canopy_shelters_the_wind) {
    // No canopy, or a sparse one (crown fill <= 5 %): the unsheltered factor from fuel depth
    CHECK(std::fabs(windAdjustmentFactor(1.0, 0, 0, 0) - unshelteredWaf(1.0)) < 1e-12);
    CHECK(std::fabs(windAdjustmentFactor(1.0, 10, 60, 0.5) - unshelteredWaf(1.0)) < 1e-12);
    // 60 % cover, crown ratio 0.5, 60 ft trees: f = 0.1, WAF = 0.555 / (sqrt(6) ln(41.6 / 7.8)) = 0.135
    const double sheltered = windAdjustmentFactor(1.0, 60, 60, 0.5);
    CHECK(std::fabs(sheltered - 0.135) < 0.002);
    CHECK(sheltered < unshelteredWaf(1.0));
}

TEST(crown_fire_thresholds) {
    // Van Wagner (1977): CBH 5 m, FMC 100 % needs ~1880 kW/m, a 2.5 m flame (RMRS-GTR-220 Table 6 / p. 24)
    const double i0 = criticalSurfaceIntensity(5.0, 100.0);
    CHECK(std::fabs(i0 - 1883.0) < 5.0);
    CHECK(std::fabs(flameLengthM(i0) - 2.5) < 0.05);
    CHECK(criticalSurfaceIntensity(10.0, 100.0) > i0);   // Higher canopy base: harder to torch
    CHECK(criticalSurfaceIntensity(5.0, 80.0) < i0);     // Drier foliage: easier
    CHECK(std::fabs(criticalActiveRos(0.15) - 20.0) < 1e-9);   // 3.0 / CBD m/min
}

namespace {
    // Uniform forest of one fuel model under the given canopy; returns how far (cells) the head fire runs
    // east in 60 minutes with a 30 km/h wind from the west, and the fire type it assigned (crown_fire.h).
    std::pair<int, int> forestRun(int fuel, float cc, float chM, float cbhM, float cbd) {
        std::vector<Cell> grid(N);
        for (Cell& c : grid) {
            c.fuelModel = static_cast<std::int16_t>(fuel);
            c.canopyCoverPct = cc; c.canopyHeightM = chM; c.canopyBaseHeightM = cbhM; c.canopyBulkDensity = cbd;
        }
        FireWeather w;
        w.windFromDeg = 270.0f;
        w.wind20ftKmh = 30.0f;
        FireSimulation sim(grid, w, 30.0f, 30.0f);
        sim.ignite(at(100, 250));
        for (int t = 0; t < 360; ++t) sim.step();      // 60 minutes
        int k = 0;
        while (sim.fireMask()[at(101 + k, 250)]) ++k;
        return { k, sim.grid()[at(100, 250)].fireType };
    }
}

TEST(crown_fire_in_the_simulation) {
    // TU5 (timber litter with shrub understory, ~640 kW/m here). A 10 m canopy base needs ~5300 kW/m, so the
    // fire stays on the ground...
    const auto [surfaceRun, surfaceType] = forestRun(165, 50.0f, 20.0f, 10.0f, 0.25f);
    CHECK(surfaceType == static_cast<int>(FireType::Surface));
    // ...a 1 m base (~170 kW/m) torches, and with dense canopy fuel (crown rate ~17 m/min > 3.0 / 0.25 = 12)
    // it becomes an active crown fire that runs far faster
    const auto [crownRun, crownType] = forestRun(165, 50.0f, 20.0f, 1.0f, 0.25f);
    CHECK(crownType == static_cast<int>(FireType::Active));
    CHECK(crownRun > 5 * std::max(surfaceRun, 1));
    // Sparse canopy fuel (CBD 0.05 needs 60 m/min) can only torch, so it spreads at the surface rate
    const auto [passiveRun, passiveType] = forestRun(165, 50.0f, 20.0f, 1.0f, 0.05f);
    CHECK(passiveType == static_cast<int>(FireType::Passive));
    CHECK(passiveRun < crownRun / 5);
    std::printf("  (60 min runs: surface %d, passive %d, active crown %d cells)\n", surfaceRun, passiveRun, crownRun);

    // Crown fire switched off: same landscape stays a surface fire
    std::vector<Cell> grid(N);
    for (Cell& c : grid) { c.fuelModel = 165; c.canopyCoverPct = 50; c.canopyHeightM = 20; c.canopyBaseHeightM = 1; c.canopyBulkDensity = 0.25f; }
    FireWeather off;
    off.crownFire = false;
    FireSimulation sim(grid, off);
    CHECK(sim.grid()[at(10, 10)].fireType == static_cast<int>(FireType::Surface));
}

TEST(canopy_slows_surface_fire) {
    // Same grass under a closed canopy (no crown fire possible: CBD 0) runs slower than in the open
    const auto [open, t1] = forestRun(102, 0.0f, 0.0f, 0.0f, 0.0f);
    const auto [shaded, t2] = forestRun(102, 80.0f, 25.0f, 5.0f, 0.0f);
    CHECK(shaded < open);
    (void)t1; (void)t2;
}

// --- Zoom and pan (viewport) ---
TEST(viewport_zooms_about_the_cursor) {
    Viewport v(1000.0f, 1000.0f);
    CHECK(v.zoom() == 1.0f && v.width() == 1000.0f);
    // Zooming in about (300, 700): that point stays at the same place on screen (30 % across, 70 % down)
    v.zoomAt(300.0f, 700.0f, 4.0f);
    CHECK(v.zoom() == 4.0f && v.width() == 250.0f);
    CHECK(std::fabs((300.0f - v.left()) / v.width() - 0.3f) < 1e-5f);
    CHECK(std::fabs((700.0f - v.top()) / v.height() - 0.7f) < 1e-5f);
    // Can't zoom out past the whole map, or in past VIEW_MAX_ZOOM
    v.zoomAt(500.0f, 500.0f, 0.01f);
    CHECK(v.zoom() == 1.0f && v.centreX() == 500.0f && v.centreY() == 500.0f);
    v.zoomAt(500.0f, 500.0f, 1000.0f);
    CHECK(v.zoom() == VIEW_MAX_ZOOM);
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
    CHECK(std::fabs(v.width() - 41.0f) < 1e-3f);
    CHECK(std::fabs(v.centreX() - 240.5f) < 1e-3f && std::fabs(v.centreY() - 250.5f) < 1e-3f);
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
    computeTopography(grid, 30.0f, 30.0f, false);
    const Cell& c = grid[at(250, 250)];
    CHECK(std::fabs(c.slopeDeg - 45.0f) < 0.01f);
    CHECK(std::fabs(c.upslopeX) < 1e-5f && std::fabs(c.upslopeY + 1.0f) < 1e-5f);
}

// --- Fire spread on the grid ---
namespace {
    FireWeather calm() {
        FireWeather w;
        w.wind20ftKmh = 0.0f;
        return w;
    }
    // Furthest cell burned along a row/column from (x0, y0), in cells
    int reach(const FireSimulation& sim, int x0, int y0, int dx, int dy) {
        int k = 0;
        while (true) {
            const int x = x0 + (k + 1) * dx, y = y0 + (k + 1) * dy;
            if (x < 0 || x >= GRID_WIDTH || y < 0 || y >= GRID_HEIGHT || !sim.fireMask()[at(x, y)]) return k;
            ++k;
        }
    }
}

// No wind, flat, uniform GR9: the front must advance at R_0 (along rows and along diagonals)
TEST(fire_front_moves_at_rothermel_rate) {
    std::vector<Cell> grid(N);
    for (Cell& c : grid) c.fuelModel = 109;
    FireSimulation sim(grid, calm(), 30.0f, 30.0f);
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
        FireSimulation sim(grid, w, 30.0f, 30.0f);
        sim.ignite(at(250, 250));
        for (int t = 0; t < 1500; ++t) sim.step();
        CHECK(reach(sim, 250, 250, 1, 0) > 3 * std::max(1, reach(sim, 250, 250, -1, 0)));
    }
    // No wind, ground rising to the north: the fire runs uphill
    {
        std::vector<Cell> grid(N);
        for (int y = 0; y < GRID_HEIGHT; ++y)
            for (int x = 0; x < GRID_WIDTH; ++x) grid[at(x, y)].elevation = (GRID_HEIGHT - y) * 15.0f;
        computeTopography(grid, 30.0f, 30.0f, false);
        FireSimulation sim(grid, calm(), 30.0f, 30.0f);
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
    FireSimulation sim(grid, w);
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
    for (int i = 0; i < GRID_WIDTH; ++i) grid[at(i, GRID_HEIGHT - 1 - i)].roadClass = 22;   // Tertiary road, SW-NE

    FireWeather w;
    w.windFromDeg = 315.0f;                          // From the NW, straight at the road
    FireSimulation sim(grid, w);
    // Positions relative to the map size, so the test holds for any GRID_WIDTH (square grid assumed):
    // k + k = GRID_WIDTH - 2 is the diagonal just north-west of the road
    const int k = (GRID_WIDTH - 2) / 2;
    CHECK(!sim.canBurn(at(k + 1, k)));               // On the road (x + y = GRID_WIDTH - 1)
    CHECK(sim.ignite(at(k - 49, k - 49)));           // North-west of the road, ~50 cells from it
    for (int t = 0; t < 4000; ++t) sim.step();
    int beyond = 0;
    for (int y = 0; y < GRID_HEIGHT; ++y)
        for (int x = 0; x < GRID_WIDTH; ++x)
            if (x + y > GRID_WIDTH - 1 && sim.fireMask()[at(x, y)]) ++beyond;
    CHECK(sim.fireMask()[at(k, k)] == 1);           // The fire reached the road...
    CHECK(beyond == 0);                              // ...and didn't cross it
}

TEST(fire_is_deterministic) {
    std::vector<Cell> grid(N);
    FireSimulation a(grid), b(grid);
    a.ignite(at(250, 250)); b.ignite(at(250, 250));
    for (int t = 0; t < 600; ++t) { a.step(); b.step(); }
    CHECK(a.fireMask() == b.fireMask());
    // The tick length only changes how often we look, not where the fire goes
    FireSimulation c(grid, FireWeather(), DEFAULT_CELL_SIZE_M, DEFAULT_CELL_SIZE_M, SIM_MINUTES_PER_TICK * 4);
    c.ignite(at(250, 250));
    for (int t = 0; t < 150; ++t) c.step();
    CHECK(c.fireMask() == a.fireMask());
}

// ============================================================================
// Planning algorithms
// ============================================================================
TEST(aco_returns_a_permutation_and_a_good_tour) {
    // Points round a circle: the optimal tour is the circle itself
    std::vector<Vec2> pts;
    const int n = 24;
    for (int i = 0; i < n; ++i) {
        const float a = 6.2831853f * i / n;
        pts.emplace_back(250.0f + 100.0f * std::cos(a), 250.0f + 100.0f * std::sin(a));
    }
    std::shuffle(pts.begin(), pts.end(), std::mt19937(5));
    std::mt19937 rng(1);
    const std::vector<int> order = solveACO(pts, rng);

    std::vector<int> sorted = order;
    std::sort(sorted.begin(), sorted.end());
    std::vector<int> expect(n);
    std::iota(expect.begin(), expect.end(), 0);
    CHECK(sorted == expect);

    float len = 0.0f;
    for (int i = 0; i < n; ++i) len += std::sqrt(sqDist(pts[order[i]], pts[order[(i + 1) % n]]));
    const float optimal = 2.0f * n * 100.0f * std::sin(3.14159265f / n);
    CHECK(len < optimal * 1.05f);
}

TEST(voronoi_assigns_every_cell) {
    std::vector<int> cells;
    for (int y = 100; y < 200; ++y) for (int x = 100; x < 200; ++x) cells.push_back(at(x, y));
    std::mt19937 rng(3);
    const Voronoi v = buildVoronoi(cells, 9, rng);
    CHECK(v.assign.size() == cells.size());
    std::vector<int> counts(9, 0);
    for (int a : v.assign) { CHECK(a >= 0 && a < 9); if (a >= 0 && a < 9) ++counts[a]; }
    // Lloyd relaxation should leave regions roughly even (each within 2x of the mean)
    for (int c : counts) CHECK(c > static_cast<int>(cells.size()) / 18 && c < static_cast<int>(cells.size()) * 2 / 9);
}

// Plans a mission around a fire with two lobes, so straight legs between sweep lanes often cross fire
TEST(planned_routes_never_cross_the_hard_no_fly_zone) {
    std::vector<std::uint8_t> fire = discFire(220, 250, 25);
    const auto lobe = discFire(290, 230, 18);
    for (int i = 0; i < N; ++i) fire[i] |= lobe[i];

    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
    live.compute(fire);
    Perimeter perim(fire);
    SafePlanner sp(live);
    std::mt19937 rng(2024);
    const MissionPlan plan = planMission(live, perim.distance(), sp, rng);

    auto checkLoop = [&](const std::vector<Vec2>& r) {
        for (size_t i = 0; i < r.size(); ++i) CHECK(sp.lineClear(r[i], r[(i + 1) % r.size()])); // incl. wrap
    };
    int nonEmpty = 0;
    for (const auto& r : plan.sweepRoutes) { checkLoop(r); nonEmpty += r.empty() ? 0 : 1; }
    checkLoop(plan.acoTour);
    CHECK(nonEmpty == SWEEP_DRONES);
    CHECK(plan.acoTour.size() >= 3);
}

// ============================================================================
// Drones
// ============================================================================
TEST(drone_moves_at_drone_speed) {
    Drone d;
    d.pos = Vec2(10.5f, 10.5f);
    d.route = { Vec2(100.5f, 10.5f), Vec2(10.5f, 10.5f) };
    stepDrone(d);
    CHECK(std::fabs(std::sqrt(sqDist(d.pos, Vec2(10.5f, 10.5f))) - DRONE_SPEED) < 1e-4f);
}

TEST(drone_stops_before_a_blocked_segment) {
    DistanceTransform edt(GRID_WIDTH, GRID_HEIGHT);
    edt.compute(discFire(60, 10, 2));     // Fire right across the drone's path
    SafePlanner sp(edt);
    Drone d;
    d.pos = Vec2(10.5f, 10.5f);
    d.route = { Vec2(100.5f, 10.5f) };
    CHECK(!stepDrone(d, &sp));            // Refuses to fly
    CHECK(sqDist(d.pos, Vec2(10.5f, 10.5f)) < 1e-6f);
}

// End-to-end: a real growing fire plus the mission controller. No drone may ever fly over a fire cell,
// and none may be inside the hard no-fly zone at the end of a frame.
TEST(mission_drones_never_fly_over_fire) {
    std::vector<Cell> grid(N);
    FireSimulation sim(grid);                       // Default weather: wind from the SW
    // Several fires so drones have to work round them
    sim.ignite(at(250, 250)); sim.ignite(at(180, 300)); sim.ignite(at(320, 200));

    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
    for (int t = 0; t < 400; ++t) sim.step();      // Let the fire get going before drones launch
    live.compute(sim.fireMask());
    Perimeter perim(sim.fireMask());
    MissionController mission(live, perim.distance(), Vec2(2.5f, GRID_HEIGHT - 2.5f));

    int crossings = 0, insideHardZone = 0;
    for (int frame = 0; frame < 1500; ++frame) {
        sim.step();
        if (sim.takeFireGrew()) {
            live.compute(sim.fireMask());
            perim.update(sim.fireMask());
            mission.markFireChanged();
        }
        std::vector<Vec2> before;
        for (const Drone& d : mission.drones()) before.push_back(d.pos);
        mission.update();
        for (size_t i = 0; i < before.size(); ++i) {
            const Vec2 p = mission.drones()[i].pos;
            if (crossesFire(sim.fireMask(), before[i], p)) ++crossings;
            const int x = static_cast<int>(p.x), y = static_cast<int>(p.y);
            if (live.squaredDistance(at(x, y)) < HARD_NO_FLY_SQ - 1e-6) ++insideHardZone;
        }
    }
    std::printf("  (drone moves over fire: %d, frames ending inside hard zone: %d)\n", crossings, insideHardZone);
    CHECK(crossings == 0);
    CHECK(insideHardZone == 0);
}

TEST(aco_drones_are_spread_round_the_tour) {
    std::vector<std::uint8_t> fire = discFire(250, 250, 30);
    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
    live.compute(fire);
    Perimeter perim(fire);
    MissionController mission(live, perim.distance(), Vec2(2.5f, GRID_HEIGHT - 2.5f));
    mission.markFireChanged();
    mission.update();   // First frame triggers the plan

    // All ACO drones start at the same base, yet must be aiming for different points on the tour
    std::set<size_t> targets;
    for (const Drone& d : mission.drones()) if (d.isACO) targets.insert(d.target);
    CHECK(static_cast<int>(targets.size()) == ACO_DRONES);
}

int main() {
    for (auto& [name, fn] : registry()) {
        const int before = g_failures;
        std::printf("%s\n", name);
        fn();
        std::printf("  %s\n", g_failures == before ? "ok" : "FAILED");
    }
    std::printf("\n%s: %d failure(s)\n", g_failures ? "FAILED" : "PASSED", g_failures);
    return g_failures ? 1 : 0;
}
