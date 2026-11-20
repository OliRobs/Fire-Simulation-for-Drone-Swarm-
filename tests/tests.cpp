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
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "aco.h"
#include "boustrophedon.h"
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

    // Centre of the grid. Tests with large fires (whose search area reaches ~2 km past the fire) put them here
    // rather than at fixed cells, so they fit on any grid size in config.h, square or not (e.g. 2000 x 2000 or
    // 1200 x 1000 cells).
    constexpr int CX = GRID_WIDTH / 2;
    constexpr int CY = GRID_HEIGHT / 2;

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

// Range-limited transform: exact within maxDistance of the sources, "very far" beyond, and correct when the
// sources move (old window cleared) or disappear altogether
TEST(edt_range_limit_matches_full_transform) {
    const int w = 300, h = 200, range = 25;
    DistanceTransform full(w, h), limited(w, h, range);
    std::mt19937 rng(11);
    for (int round = 0; round < 3; ++round) {
        std::vector<std::uint8_t> src(w * h, 0);
        // Round 0: sources bottom-right; round 1: top-left only (window moves); round 2: no sources
        const int ox = round == 0 ? 180 : 20, oy = round == 0 ? 120 : 10;
        if (round < 2) for (int k = 0; k < 8; ++k) src[(oy + rng() % 40) * w + ox + rng() % 60] = 1;
        full.compute(src);
        limited.compute(src);
        for (int i = 0; i < w * h; ++i) {
            const double f = full.squaredDistance(i), l = limited.squaredDistance(i);
            if (f <= double(range) * range) CHECK(std::fabs(f - l) < 1e-9);   // Exact inside the range
            else CHECK(l >= double(range) * range);                          // Never nearer than it really is
        }
    }
}

// Giving compute() a rectangle that holds the sources (as main.cpp does with the fire's bounding box) must give
// the same answer as letting it scan for them, including when the rectangle is looser than the sources
TEST(edt_with_source_bounds_matches_scan) {
    const int w = 300, h = 200, range = 25;
    DistanceTransform scanned(w, h, range), hinted(w, h, range), full(w, h);
    std::mt19937 rng(12);
    for (int round = 0; round < 3; ++round) {
        std::vector<std::uint8_t> src(w * h, 0);
        GridRect box;   // Exact bounding box of the sources
        const int ox = round == 0 ? 180 : 20, oy = round == 0 ? 120 : 10;
        if (round < 2) for (int k = 0; k < 8; ++k) {
            const int x = ox + rng() % 60, y = oy + rng() % 40;
            src[y * w + x] = 1;
            box.include(x, y);
        }
        // Loosen the box by a few cells (and let it run off the grid on round 1) to stand in for a bigger fire
        const GridRect loose{ box.x0 - 30, box.y0 - 15, box.x1 + 7, box.y1 + 3 };
        scanned.compute(src);
        hinted.compute(src, round == 2 ? GridRect{} : loose);   // Round 2: no sources, empty rectangle
        full.compute(src);
        const GridRect win = hinted.window();
        for (int i = 0; i < w * h; ++i) {
            const double f = full.squaredDistance(i), s = scanned.squaredDistance(i), l = hinted.squaredDistance(i);
            if (f <= double(range) * range) {
                CHECK(std::fabs(f - l) < 1e-9 && std::fabs(s - l) < 1e-9);   // Exact inside the range
                const int x = i % w, y = i / w;
                CHECK(x >= win.x0 && x < win.x1 && y >= win.y0 && y < win.y1); // ...and inside window()
            }
            else CHECK(l >= double(range) * range);                          // Never nearer than it really is
        }
    }
}

// The fire simulation's bounding box must always be exactly the box of its fire mask
TEST(fire_bounds_track_the_fire_mask) {
    std::vector<Cell> grid(N);
    FireSimulation sim(grid);
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

// Perimeter::update with the fire's rectangle (main.cpp) and without it (scanning) must move the perimeter
// identically: same ring, same distances everywhere the mission planner reads them
TEST(perimeter_update_with_bounds_matches_scan) {
    std::vector<Cell> grid(N);
    FireSimulation sim(grid);
    sim.ignite(at(400, 400));
    for (int t = 0; t < 200; ++t) sim.step();
    Perimeter withBounds(sim.fireMask()), scanned(sim.fireMask());
    int moves = 0;
    for (int t = 0; t < 1500; ++t) {
        sim.step();
        if (!sim.takeFireGrew()) continue;
        const bool a = withBounds.update(sim.fireMask(), sim.fireBounds());
        const bool b = scanned.update(sim.fireMask());
        CHECK(a == b);
        moves += a ? 1 : 0;
    }
    CHECK(moves > 0);   // The test must actually exercise the perimeter moving
    CHECK(withBounds.ringCells() == scanned.ringCells());
    constexpr double reach = (SEARCH_DISTANCE + 0.5) * (SEARCH_DISTANCE + 0.5);
    for (int i = 0; i < N; ++i) {
        const double a = withBounds.distance().squaredDistance(i), b = scanned.distance().squaredDistance(i);
        if (a < reach || b < reach) CHECK(std::fabs(a - b) < 1e-9);
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
    // The view keeps the map's aspect ratio, so on a non-square map only the tighter side is exactly 41 cells;
    // the other is wider, and the whole 41 x 41 box is always in view
    CHECK(std::fabs(std::min(v.width(), v.height()) - 41.0f) < 1e-3f);
    CHECK(v.width() >= 41.0f - 1e-3f && v.height() >= 41.0f - 1e-3f);
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
    computeTopography(coarse, 10.0f, 10.0f, false);
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
    // The road runs corner to corner across the largest square in the top-left of the map, from (0, S - 1) to
    // (S - 1, 0) with S the shorter side, so it meets two map edges and cuts the map in two whatever its shape
    const int S = std::min(GRID_WIDTH, GRID_HEIGHT);
    for (int i = 0; i < S; ++i) grid[at(i, S - 1 - i)].roadClass = 22;   // Tertiary road, SW-NE

    FireWeather w;
    w.windFromDeg = 315.0f;                          // From the NW, straight at the road
    FireSimulation sim(grid, w);
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

// solveACO builds each ant's tour from a per-iteration choice table and a list of unvisited points. This
// checks it against the original version (all n points scanned per step, tau * heur worked out each time),
// kept here unchanged: the tours and the random numbers drawn must be identical, bit for bit.
TEST(aco_choice_table_matches_original_scan) {
    // --- Reference: the original solveACO, copied unchanged (only renamed) ---
    auto reference = [](const std::vector<Vec2>& pts, std::mt19937& rng) {
        const int n = static_cast<int>(pts.size());
        std::vector<int> best(n);
        std::iota(best.begin(), best.end(), 0);
        if (n <= 3) return best;
        std::vector<float> d(static_cast<size_t>(n) * n), heur(static_cast<size_t>(n) * n, 0.0f);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) {
                d[i * n + j] = std::sqrt(sqDist(pts[i], pts[j]));
                if (i != j) { const float h = 1.0f / (d[i * n + j] + 1e-3f); heur[i * n + j] = h * h * h; }
            }
        auto tourLength = [&](const std::vector<int>& t) {
            float len = 0.0f;
            for (int i = 0; i < n; ++i) len += d[t[i] * n + t[(i + 1) % n]];
            return len;
        };
        {
            std::vector<char> seen(n, 0);
            int cur = 0;
            best[0] = 0; seen[0] = 1;
            for (int step = 1; step < n; ++step) {
                int next = -1; float nd = std::numeric_limits<float>::max();
                for (int j = 0; j < n; ++j) if (!seen[j] && d[cur * n + j] < nd) { nd = d[cur * n + j]; next = j; }
                best[step] = next; seen[next] = 1; cur = next;
            }
        }
        float bestLen = tourLength(best);
        std::vector<float> tau(static_cast<size_t>(n) * n, ACO_ANTS / std::max(bestLen, 1e-3f));
        std::vector<std::vector<int>> tours(ACO_ANTS, std::vector<int>(n));
        std::vector<float> lens(ACO_ANTS), w(n);
        std::vector<char> visited(n);
        std::uniform_real_distribution<float> u01(0.0f, 1.0f);
        for (int it = 0; it < ACO_ITERATIONS; ++it) {
            for (int a = 0; a < ACO_ANTS; ++a) {
                std::fill(visited.begin(), visited.end(), 0);
                int cur = static_cast<int>(rng() % n);
                tours[a][0] = cur; visited[cur] = 1;
                for (int step = 1; step < n; ++step) {
                    float sum = 0.0f;
                    for (int j = 0; j < n; ++j) {
                        w[j] = visited[j] ? 0.0f : tau[cur * n + j] * heur[cur * n + j];
                        sum += w[j];
                    }
                    int next = -1;
                    if (sum > 0.0f) {
                        float pick = u01(rng) * sum;
                        for (int j = 0; j < n; ++j) {
                            if (w[j] <= 0.0f) continue;
                            next = j; pick -= w[j];
                            if (pick <= 0.0f) break;
                        }
                    }
                    if (next < 0) for (int j = 0; j < n; ++j) if (!visited[j]) { next = j; break; }
                    tours[a][step] = next; visited[next] = 1; cur = next;
                }
                lens[a] = tourLength(tours[a]);
                if (lens[a] < bestLen) { bestLen = lens[a]; best = tours[a]; }
            }
            for (float& t : tau) t *= (1.0f - ACO_EVAPORATION);
            auto deposit = [&](const std::vector<int>& t, float amount) {
                for (int i = 0; i < n; ++i) {
                    const int a = t[i], b = t[(i + 1) % n];
                    tau[a * n + b] += amount;
                    tau[b * n + a] += amount;
                }
            };
            for (int a = 0; a < ACO_ANTS; ++a) deposit(tours[a], 1.0f / lens[a]);
            deposit(best, ACO_ELITE_WEIGHT / bestLen);
        }
        return best;
    };

    // --- Point sets: random spreads of many sizes, clusters with repeated points, and points so far apart
    // that the weights underflow to zero (this exercises the deterministic fallback pick) ---
    std::mt19937 gen(11);
    std::uniform_real_distribution<float> coord(0.0f, 2000.0f);
    std::vector<std::vector<Vec2>> sets;
    for (int n : { 4, 5, 9, 24, 60, 130 }) {
        std::vector<Vec2> pts;
        for (int i = 0; i < n; ++i) pts.emplace_back(coord(gen), coord(gen));
        sets.push_back(pts);
    }
    {
        std::vector<Vec2> pts;   // Clusters of points, some exactly on top of each other
        for (int i = 0; i < 40; ++i) pts.emplace_back(100.0f * (i % 4), 100.0f * ((i / 4) % 3));
        sets.push_back(pts);
    }
    {
        std::vector<Vec2> pts;   // 1e15 cells apart: 1 / d^3 underflows, every weight is 0
        for (int i = 0; i < 12; ++i) pts.emplace_back(1e15f * i, 1e15f * (i % 3));
        sets.push_back(pts);
    }

    // --- Same seed into both: same tour, and the random stream left in the same place ---
    int mismatches = 0;
    for (size_t s = 0; s < sets.size(); ++s) {
        std::mt19937 r1(100 + static_cast<unsigned>(s)), r2(100 + static_cast<unsigned>(s));
        const std::vector<int> got = solveACO(sets[s], r1);
        const std::vector<int> want = reference(sets[s], r2);
        if (got != want || r1() != r2()) ++mismatches;
    }
    CHECK(mismatches == 0);
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

// The bucketed nearest-site search inside buildVoronoi must give exactly what testing every site gives
// (lowest site index on a tie), on an awkward shape: a ring with a gap plus scattered islands of cells
TEST(voronoi_bucketed_assignment_matches_brute_force) {
    std::vector<int> cells;
    std::mt19937 shape(8);
    for (int y = 300; y < 700; ++y)
        for (int x = 300; x < 700; ++x) {
            const int r2 = (x - 500) * (x - 500) + (y - 500) * (y - 500);
            const bool ring = r2 > 120 * 120 && r2 < 190 * 190 && !(x > 480 && x < 530 && y < 500);
            if (ring || shape() % 97 == 0) cells.push_back(at(x, y));
        }
    for (int sites : { 1, 7, 64 }) {
        std::mt19937 rng(4);
        const Voronoi v = buildVoronoi(cells, sites, rng);
        CHECK(static_cast<int>(v.sites.size()) == sites);
        int wrong = 0;
        for (size_t i = 0; i < cells.size(); ++i) {
            const Vec2 p = cellCentre(cells[i]);
            int best = 0;
            float bd = std::numeric_limits<float>::max();
            for (int k = 0; k < sites; ++k) {
                const float dd = sqDist(p, v.sites[k]);
                if (dd < bd) { bd = dd; best = k; }
            }
            wrong += v.assign[i] != best ? 1 : 0;
        }
        CHECK(wrong == 0);
    }
}

// The bucketed farthest-point seeding inside buildVoronoi must pick exactly the seeds the plain method picks
// (update every cell's nearest-seed distance per seed, take the first farthest cell). The reference below is
// the whole algorithm done the plain way: plain seeding, then Lloyd passes testing every site for every cell.
// Sites and labels must match bit for bit, on shapes with many equal distances (a solid square, where ties
// between cells are common) and with gaps (a ring with a gap plus scattered islands).
TEST(voronoi_bucketed_seeding_matches_plain_method) {
    // Plain reference implementation of buildVoronoi, written as simply as possible
    auto plainVoronoi = [](const std::vector<int>& cells, int S, std::mt19937& rng) {
        const size_t n = cells.size();
        std::vector<Vec2> pts(n), sites;
        for (size_t i = 0; i < n; ++i) pts[i] = cellCentre(cells[i]);
        // Farthest-point seeding: first seed random, then the first cell farthest from all seeds so far
        std::vector<float> minD(n, std::numeric_limits<float>::max());
        sites.push_back(pts[rng() % n]);
        while (static_cast<int>(sites.size()) < S) {
            size_t far = 0;
            float farD = -1.0f;
            for (size_t i = 0; i < n; ++i) {
                minD[i] = std::min(minD[i], sqDist(pts[i], sites.back()));
                if (minD[i] > farD) { farD = minD[i]; far = i; }
            }
            sites.push_back(pts[far]);
        }
        // Nearest site for every cell (lowest index on a tie), then Lloyd moves each site to its mean
        std::vector<int> assign(n, 0);
        auto assignAll = [&]() {
            for (size_t i = 0; i < n; ++i) {
                float bd = std::numeric_limits<float>::max();
                for (int k = 0; k < S; ++k) {
                    const float dd = sqDist(pts[i], sites[k]);
                    if (dd < bd) { bd = dd; assign[i] = k; }
                }
            }
        };
        for (int it = 0; it < LLOYD_ITERATIONS; ++it) {
            assignAll();
            std::vector<Vec2> sum(S, Vec2(0.0f, 0.0f));
            std::vector<int> cnt(S, 0);
            for (size_t i = 0; i < n; ++i) { sum[assign[i]] += pts[i]; ++cnt[assign[i]]; }
            for (int k = 0; k < S; ++k) if (cnt[k] > 0) sites[k] = sum[k] / static_cast<float>(cnt[k]);
        }
        assignAll();
        Voronoi v;
        v.sites = sites;
        v.assign = assign;
        return v;
    };

    // Shape 1: solid square (lots of exactly equal distances, so the tie-break matters)
    std::vector<int> square;
    for (int y = 400; y < 520; ++y) for (int x = 600; x < 720; ++x) square.push_back(at(x, y));
    // Shape 2: ring with a gap plus scattered islands, spread over many mostly empty blocks
    std::vector<int> ring;
    std::mt19937 shape(11);
    for (int y = 200; y < 800; ++y)
        for (int x = 200; x < 800; ++x) {
            const int r2 = (x - 500) * (x - 500) + (y - 500) * (y - 500);
            const bool onRing = r2 > 200 * 200 && r2 < 280 * 280 && !(y > 480 && y < 540 && x > 500);
            if (onRing || shape() % 131 == 0) ring.push_back(at(x, y));
        }

    for (const auto* cells : { &square, &ring }) {
        for (int sites : { 1, 2, 9, 64 }) {
            std::mt19937 rngA(77), rngB(77);
            const Voronoi fast = buildVoronoi(*cells, sites, rngA);
            const Voronoi ref = plainVoronoi(*cells, sites, rngB);
            CHECK(fast.sites.size() == ref.sites.size());
            bool sameSites = fast.sites.size() == ref.sites.size();
            for (size_t k = 0; sameSites && k < ref.sites.size(); ++k)
                sameSites = fast.sites[k].x == ref.sites[k].x && fast.sites[k].y == ref.sites[k].y;
            CHECK(sameSites);
            CHECK(fast.assign == ref.assign);
        }
    }
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
    CHECK(nonEmpty == static_cast<int>(plan.sweepRoutes.size()));   // Every sector has a route
    CHECK(nonEmpty >= SWEEP_DRONES);
    CHECK(plan.acoTour.size() >= 3);
}

// ============================================================================
// Perimeter line (ordered loops traced from the committed-fire distance field)
// ============================================================================
TEST(perimeter_line_round_fire_is_one_loop) {
    const int r = 30;
    Perimeter perim(discFire(400, 400, r));
    const auto loops = tracePerimeterLine(perim.distance());
    CHECK(loops.size() == 1);
    if (loops.size() != 1) return;
    // About 2 pi (r + 1.2 km) long, and every point on it 1.2 km from the fire's edge
    const float expected = 2.0f * 3.14159265f * (r + PERIMETER_DISTANCE);
    CHECK(std::fabs(loopLength(loops[0]) - expected) < 0.02f * expected);
    float worst = 0.0f;
    for (const Vec2& p : loops[0]) {
        const float d = std::sqrt(sqDist(p, Vec2(400.5f, 400.5f))) - (r + PERIMETER_DISTANCE);
        worst = std::max(worst, std::fabs(d));
    }
    CHECK(worst < 1.0f);
    // Anticlockwise on screen = negative signed area in grid coordinates (y down)
    CHECK(loopSignedArea(loops[0]) < 0.0);
    // Consecutive points are close together (one square apart at most), i.e. the loop really is ordered
    bool ordered = true;
    for (size_t i = 0; i < loops[0].size(); ++i)
        ordered = ordered && sqDist(loops[0][i], loops[0][(i + 1) % loops[0].size()]) < 2.01f;
    CHECK(ordered);
}

TEST(perimeter_line_two_fires_give_two_loops_the_same_way_round) {
    std::vector<std::uint8_t> fire = discFire(400, 400, 20);
    const auto other = discFire(GRID_WIDTH * 3 / 4, GRID_HEIGHT * 3 / 4, 35);   // Far enough apart on any grid
    for (int i = 0; i < N; ++i) fire[i] |= other[i];
    Perimeter perim(fire);
    const auto loops = tracePerimeterLine(perim.distance());
    CHECK(loops.size() == 2);
    for (const auto& l : loops) CHECK(loopSignedArea(l) < 0.0);
}

// ============================================================================
// Sweep lanes and sectors
// ============================================================================
namespace {
    // Everything a sweep test needs for one fire: the distance fields, the planner and the plan
    struct SweepSetup {
        DistanceTransform live{ GRID_WIDTH, GRID_HEIGHT };
        Perimeter perim;
        SafePlanner sp;
        MissionPlan plan;
        explicit SweepSetup(const std::vector<std::uint8_t>& fire) : perim(fire), sp(live) {
            live.compute(fire);
            std::mt19937 rng(2024);
            plan = planMission(live, perim.distance(), sp, rng);
        }
    };

    // Distance from p to the nearest run of any lane
    float nearestLane(const std::vector<SweepLane>& lanes, Vec2 p) {
        float best = std::numeric_limits<float>::max();
        for (const SweepLane& l : lanes)
            for (size_t i = 0; i + 1 < l.pts.size(); i += 2) {
                const Vec2 a = l.pts[i], ab = l.pts[i + 1] - a, ap = p - a;
                const float l2 = ab.x * ab.x + ab.y * ab.y;
                const float t = l2 > 1e-9f ? std::clamp((ap.x * ab.x + ap.y * ab.y) / l2, 0.0f, 1.0f) : 0.0f;
                best = std::min(best, std::sqrt(sqDist(p, a + ab * t)));
            }
        return best;
    }

    // Count of inner-band cells that no lane sees (further than LANE_COVER_RADIUS from every lane)
    int unseenInnerCells(const MissionPlan& plan) {
        int unseen = 0;
        for (int y = plan.bounds.y0; y < plan.bounds.y1; ++y)
            for (int x = plan.bounds.x0; x < plan.bounds.x1; ++x)
                if (plan.zone[at(x, y)] == 1 && nearestLane(plan.lanes, Vec2(x + 0.5f, y + 0.5f)) > LANE_COVER_RADIUS + 1e-3f)
                    ++unseen;
        return unseen;
    }

    // A fire with a concave bay: a disc with a bite taken out of its east side
    std::vector<std::uint8_t> bayFire() {
        std::vector<std::uint8_t> fire = discFire(500, 500, 90);
        const auto bite = discFire(575, 500, 60);
        for (int i = 0; i < N; ++i) if (bite[i]) fire[i] = 0;
        return fire;
    }
}

TEST(sweep_lanes_point_at_a_round_fire) {
    const int r = 30;
    const Vec2 c(400.5f, 400.5f);
    SweepSetup s(discFire(400, 400, r));
    int mainLanes = 0, badDir = 0, badStart = 0, badEnd = 0;
    for (const SweepLane& l : s.plan.lanes) {
        if (l.gapFill) continue;
        ++mainLanes;
        // Direction: from the lane's outer end straight at the fire's centre, to within 3 degrees
        const Vec2 along = l.pts.back() - l.pts.front(), toC = c - l.pts.front();
        const float cosA = (along.x * toC.x + along.y * toC.y)
            / (std::sqrt(along.x * along.x + along.y * along.y) * std::sqrt(toC.x * toC.x + toC.y * toC.y));
        if (cosA < std::cos(3.0f * 3.14159265f / 180.0f)) ++badDir;
        // Runs from the perimeter line (band's outer edge) to the 300 m standoff
        const float rOut = std::sqrt(sqDist(l.pts.front(), c)), rIn = std::sqrt(sqDist(l.pts.back(), c));
        if (std::fabs(rOut - (r + PERIMETER_DISTANCE)) > 2.0f) ++badStart;
        if (std::fabs(rIn - (r + MIN_FIRE_DISTANCE)) > 2.0f) ++badEnd;
    }
    std::printf("  (%d main lanes; off direction %d, off start %d, off end %d)\n", mainLanes, badDir, badStart, badEnd);
    // One lane per 120 m (at most) of perimeter line
    CHECK(mainLanes == static_cast<int>(std::ceil(loopLength(s.plan.perimeterLoops[0]) / LANE_SPACING)));
    CHECK(badDir == 0);
    CHECK(badStart == 0);
    CHECK(badEnd == 0);
}

TEST(sweep_lanes_cover_the_whole_band) {
    SweepSetup round(discFire(400, 400, 30));
    int gapLanes = 0;
    for (const SweepLane& l : round.plan.lanes) gapLanes += l.gapFill;
    CHECK(gapLanes == 0);                      // Convex fire: the main lanes alone see everything
    CHECK(unseenInnerCells(round.plan) == 0);

    SweepSetup bay(bayFire());
    gapLanes = 0;
    for (const SweepLane& l : bay.plan.lanes) gapLanes += l.gapFill;
    std::printf("  (bay: %zu lanes, %d of them gap fill)\n", bay.plan.lanes.size(), gapLanes);
    CHECK(gapLanes > 0);                       // Lanes spread apart in the bay, so gap fill had work to do
    CHECK(unseenInnerCells(bay.plan) == 0);
}

TEST(sweep_sectors_are_exclusive_and_contiguous) {
    SweepSetup s(bayFire());
    const MissionPlan& plan = s.plan;
    const int sectors = static_cast<int>(plan.sweepRoutes.size());

    // Every inner cell has exactly one owner (a valid sector); no cell outside the band has one
    int unowned = 0, strayOwner = 0;
    for (int i = 0; i < N; ++i) {
        if (plan.zone[i] == 1) unowned += plan.owner[i] < 0 || plan.owner[i] >= sectors;
        else strayOwner += plan.owner[i] >= 0;
    }
    CHECK(unowned == 0);
    CHECK(strayOwner == 0);

    // Each sector is one contiguous run of the lane sequence
    std::vector<int> runsOf(sectors, 0);
    for (size_t i = 0; i < plan.lanes.size(); ++i) {
        CHECK(plan.lanes[i].sector >= 0 && plan.lanes[i].sector < sectors);
        if (i == 0 || plan.lanes[i].sector != plan.lanes[i - 1].sector) ++runsOf[plan.lanes[i].sector];
    }
    for (int r : runsOf) CHECK(r == 1);

    // At the perimeter line, the edge lanes of neighbouring sectors are one lane spacing apart (+- 1 cell),
    // i.e. the border is halfway between them and each edge lane is 60 m inside it. Measured on a round fire,
    // where the straight distance between seeds equals the distance along the line (the bay puts a sharp
    // corner in the perimeter line, where the straight distance is shorter).
    SweepSetup round(discFire(400, 400, 30));
    std::vector<const SweepLane*> mainLanes;
    for (const SweepLane& l : round.plan.lanes) if (!l.gapFill) mainLanes.push_back(&l);
    int borders = 0;
    for (size_t i = 0; i < mainLanes.size(); ++i) {
        const SweepLane& a = *mainLanes[i];
        const SweepLane& b = *mainLanes[(i + 1) % mainLanes.size()];
        if (a.sector == b.sector) continue;
        ++borders;
        CHECK(std::fabs(std::sqrt(sqDist(a.seed, b.seed)) - LANE_SPACING) <= 1.0f);
    }
    CHECK(borders == static_cast<int>(round.plan.sweepRoutes.size()));   // One loop: a border per sector
}

// ============================================================================
// Sweep fleet sizing
// ============================================================================
// The incremental lap builder used by the sector split gives exactly (bit for bit) the lap the plain
// lane-by-lane sum gives, for every sector of a random set of lanes, whether built with addNext() or
// asked ahead with lapWithNext(). Exactness matters: the sector split compares these laps with a limit,
// so any rounding difference could move a sector border and change the plan.
TEST(sector_lap_builder_matches_plain_sum) {
    // Reference: the original estimateSectorLength loop, kept here unchanged
    auto plain = [](const std::vector<SweepLane>& lanes, int first, int last) {
        float len = 0.0f;
        Vec2 start, end;
        for (int i = first; i <= last; ++i) {
            const std::vector<Vec2>& p = lanes[i].pts;
            const bool inward = (i - first) % 2 == 0;
            const Vec2 a = inward ? p.front() : p.back(), b = inward ? p.back() : p.front();
            if (i == first) start = a;
            else len += std::sqrt(sqDist(end, a));
            len += lanes[i].length();
            end = b;
        }
        return len + std::sqrt(sqDist(end, start));
    };

    // Random lanes of 1 to 3 runs (2 to 6 points, single-cell runs included), scattered over the map
    std::mt19937 rng(7);
    std::uniform_real_distribution<float> coord(0.0f, 1000.0f);
    std::vector<SweepLane> lanes(40);
    for (SweepLane& l : lanes) {
        const int pts = 2 * (1 + static_cast<int>(rng() % 3));
        for (int k = 0; k < pts; ++k) l.pts.push_back(Vec2(coord(rng), coord(rng)));
        if (rng() % 5 == 0) l.pts[1] = l.pts[0];                     // A single-cell run
    }

    int mismatches = 0;
    for (int first = 0; first < static_cast<int>(lanes.size()); ++first) {
        SectorLapBuilder b(lanes, first);
        for (int last = first;; ++last) {
            const float want = plain(lanes, first, last);
            if (b.lap() != want || estimateSectorLength(lanes, first, last) != want) ++mismatches;
            if (last + 1 >= static_cast<int>(lanes.size())) break;
            const float ahead = b.lapWithNext();
            if (ahead != plain(lanes, first, last + 1)) ++mismatches;
            b.addNext();
            if (b.last() != last + 1 || b.lap() != ahead) ++mismatches;
        }
    }
    CHECK(mismatches == 0);
}

TEST(sweep_laps_stay_within_the_cap) {
    for (int r : { 1, 50, 100 }) {   // New fire, 1 km and 2 km across
        SweepSetup s(discFire(CX, CY, r));
        const int sectors = static_cast<int>(s.plan.sweepRoutes.size());
        std::printf("  (fire %d m across: %d sweep drones, laps", 2 * r * static_cast<int>(SIM_CELL_SIZE_M), sectors);
        for (float lap : s.plan.sweepLapS) std::printf(" %.1f", lap / 60.0f);
        std::printf(" min)\n");
        CHECK(sectors >= MIN_SWEEP_DRONES && sectors <= MAX_SWEEP_DRONES);
        CHECK(sectors < MAX_SWEEP_DRONES);   // These fires are small enough not to need the maximum
        CHECK(s.plan.sweepNeeded == sectors);
        CHECK(!s.plan.sweepOverLimit);
        for (int k = 0; k < sectors; ++k) {
            // The lap includes the closing leg back to the start
            CHECK(std::fabs(s.plan.sweepLapS[k] - flightSeconds(closedRouteLength(s.plan.sweepRoutes[k]))) < 1e-3f);
            CHECK(s.plan.sweepLapS[k] <= MAX_SWEEP_LAP_S);
        }
    }
    // Bigger fires need more sweep drones
    SweepSetup small(discFire(CX, CY, 1)), big(discFire(CX, CY, 100));
    CHECK(big.plan.sweepRoutes.size() > small.plan.sweepRoutes.size());
}

// Two re-plans with a small change to the fire keep most cells with the same drone
TEST(sweep_ownership_is_stable) {
    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
    std::vector<std::uint8_t> fire = discFire(600, 600, 60);
    live.compute(fire);
    Perimeter perim(fire);
    SafePlanner sp(live);
    std::mt19937 rng(2024);
    MissionPlan first = planMission(live, perim.distance(), sp, rng);
    // As the controller does: new drones (sectorDrone -1) get fresh ids
    int nextId = 100;
    for (int& id : first.sectorDrone) if (id < 0) id = nextId++;

    // The fire grows a little on one side
    const auto lobe = discFire(665, 600, 12);
    for (int i = 0; i < N; ++i) fire[i] |= lobe[i];
    live.compute(fire);
    perim.update(fire);
    SweepFleet fleet;
    for (size_t k = 0; k < first.sectorDrone.size(); ++k) {
        fleet.ids.push_back(first.sectorDrone[k]);
        fleet.positions.push_back(first.sweepRoutes[k].empty() ? Vec2(0.5f, 0.5f) : first.sweepRoutes[k].front());
    }
    fleet.previous = &first;
    rng.seed(2024);
    const MissionPlan second = planMission(live, perim.distance(), sp, rng, fleet);

    long long both = 0, same = 0;
    for (int i = 0; i < N; ++i) {
        if (first.owner[i] < 0 || second.owner[i] < 0) continue;
        ++both;
        same += first.sectorDrone[first.owner[i]] == second.sectorDrone[second.owner[i]];
    }
    std::printf("  (%.1f %% of %lld cells kept their drone)\n", 100.0 * same / std::max(1LL, both), both);
    CHECK(both > 0);
    CHECK(same >= 0.8 * both);
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

    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT, LIVE_EDT_RANGE);   // Range-limited, as in main.cpp
    for (int t = 0; t < 400; ++t) sim.step();      // Let the fire get going before drones launch
    live.compute(sim.fireMask());
    Perimeter perim(sim.fireMask());
    MissionController mission(live, perim.distance(), Vec2(2.5f, GRID_HEIGHT - 2.5f));

    int crossings = 0, insideHardZone = 0;
    for (int frame = 0; frame < 1500; ++frame) {
        sim.step();
        if (sim.takeFireGrew()) {
            live.compute(sim.fireMask(), sim.fireBounds());   // As main.cpp: the fire's box saves a grid scan
            perim.update(sim.fireMask(), sim.fireBounds());
            mission.markFireChanged();
        }
        // Positions by drone id: the fleet grows and shrinks, so list positions don't stay with one drone
        std::vector<std::pair<int, Vec2>> before;
        for (const Drone& d : mission.drones()) before.emplace_back(d.id, d.pos);
        mission.update();
        for (const Drone& d : mission.drones()) {
            for (const auto& [id, was] : before)
                if (id == d.id && crossesFire(sim.fireMask(), was, d.pos)) ++crossings;
            const int x = static_cast<int>(d.pos.x), y = static_cast<int>(d.pos.y);
            if (live.squaredDistance(at(x, y)) < HARD_NO_FLY_SQ - 1e-6) ++insideHardZone;
        }
    }
    std::printf("  (drone moves over fire: %d, frames ending inside hard zone: %d)\n", crossings, insideHardZone);
    CHECK(crossings == 0);
    CHECK(insideHardZone == 0);
}

// Growing a fire makes the controller launch more drones, at once, from the launch point
TEST(fleet_grows_with_the_fire) {
    const Vec2 base(2.5f, GRID_HEIGHT - 2.5f);
    std::vector<std::uint8_t> fire = discFire(CX, CY, 1);
    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
    live.compute(fire);
    Perimeter perim(fire);
    MissionController mission(live, perim.distance(), base);
    mission.markFireChanged();
    mission.update();
    const int sweepBefore = mission.activeCount(false);
    int maxIdBefore = -1;
    for (const Drone& d : mission.drones()) maxIdBefore = std::max(maxIdBefore, d.id);

    fire = discFire(CX, CY, 100);   // Now 2 km across
    live.compute(fire);
    perim.update(fire);
    mission.markFireChanged();
    for (int f = 0; f < REPLAN_INTERVAL; ++f) mission.update();
    const int sweepAfter = mission.activeCount(false);
    std::printf("  (sweep drones %d -> %d)\n", sweepBefore, sweepAfter);
    CHECK(sweepAfter > sweepBefore);
    CHECK(sweepAfter == mission.plan().sweepNeeded);
    // Every sector has exactly one active drone flying it
    std::set<int> flying;
    for (const Drone& d : mission.drones()) if (!d.isACO && d.state == DroneState::Active) flying.insert(d.id);
    CHECK(static_cast<int>(flying.size()) == static_cast<int>(mission.plan().sweepRoutes.size()));
    for (int id : mission.plan().sectorDrone) CHECK(flying.count(id) == 1);
    // The new drones have new ids
    int newOnes = 0;
    for (const Drone& d : mission.drones()) newOnes += d.id > maxIdBefore;
    CHECK(newOnes == sweepAfter - sweepBefore);
}

// A shrinking fire lowers each count by one only after SHRINK_REPLANS re-plans in a row; retired drones fly
// home without ever crossing the fire and are removed there. Then a growing fire raises both counts at once.
TEST(fleet_shrinks_slowly_and_retired_drones_fly_home) {
    const Vec2 base(2.5f, GRID_HEIGHT - 2.5f);
    std::vector<std::uint8_t> fire = discFire(CX, CY, 100);   // 2 km across
    // Plain distance fields (rather than a Perimeter, whose committed fire never shrinks) so the test can
    // make the fire smaller. The committed-fire field has the Perimeter's range limit, as the planner expects.
    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT), perimDt(GRID_WIDTH, GRID_HEIGHT, SEARCH_DISTANCE + 1);
    live.compute(fire);
    perimDt.compute(fire);
    MissionController mission(live, perimDt, base);

    int crossings = 0, insideHardZone = 0;
    std::set<int> retired;          // Every drone ever seen returning
    auto frames = [&](int count) {
        for (int f = 0; f < count; ++f) {
            std::vector<std::pair<int, Vec2>> before;
            for (const Drone& d : mission.drones()) before.emplace_back(d.id, d.pos);
            mission.update();
            for (const Drone& d : mission.drones()) {
                if (d.state == DroneState::Returning) retired.insert(d.id);
                for (const auto& [id, was] : before)
                    if (id == d.id && crossesFire(fire, was, d.pos)) ++crossings;
                if (live.squaredDistance(at(static_cast<int>(d.pos.x), static_cast<int>(d.pos.y))) < HARD_NO_FLY_SQ - 1e-6)
                    ++insideHardZone;
            }
        }
    };
    // One re-plan: mark the plan dirty and run until the controller has re-planned
    auto replan = [&]() { mission.markFireChanged(); frames(REPLAN_INTERVAL); };
    auto setFire = [&](int r) { fire = discFire(CX, CY, r); live.compute(fire); perimDt.compute(fire); };

    replan();
    const int sweepBig = mission.activeCount(false);
    CHECK(sweepBig == 8);
    // Bring the ACO fleet down to what the big fire needs, so the counts below start from a settled fleet
    for (int k = 0; k < 3 * SHRINK_REPLANS && mission.activeCount(true) > mission.plan().acoNeeded; ++k) replan();
    const int acoBig = mission.activeCount(true);
    CHECK(acoBig == mission.plan().acoNeeded);

    // --- Shrink the fire to a new fire: each count falls by one only on the SHRINK_REPLANS-th re-plan ---
    setFire(1);
    for (int k = 1; k <= SHRINK_REPLANS; ++k) {
        replan();
        const int expectSweep = k < SHRINK_REPLANS ? sweepBig : sweepBig - 1;
        CHECK(mission.activeCount(false) == expectSweep);
        CHECK(mission.activeCount(true) == (k < SHRINK_REPLANS ? acoBig : std::max(MIN_ACO_DRONES, acoBig - 1)));
    }
    std::printf("  (after %d re-plans on a small fire: sweep %d -> %d, ACO %d -> %d)\n", SHRINK_REPLANS, sweepBig,
                mission.activeCount(false), acoBig, mission.activeCount(true));
    CHECK(!retired.empty());

    // Keep re-planning: both fleets settle. The sweep fleet stops one above the bare minimum needed when one
    // fewer drone would push a lap past SHRINK_MARGIN of the cap (that margin is the point of the rule).
    for (int k = 0; k < 4 * SHRINK_REPLANS; ++k) replan();
    CHECK(mission.activeCount(false) >= mission.plan().sweepNeeded);
    CHECK(mission.activeCount(false) == mission.plan().sweepNeeded || !mission.plan().sweepShrinkOk);
    for (float lap : mission.plan().sweepLapS) CHECK(lap <= MAX_SWEEP_LAP_S);
    CHECK(mission.activeCount(true) == mission.plan().acoNeeded);
    CHECK(mission.activeCount(true) == MIN_ACO_DRONES);

    // Every retired drone gets home and is removed (allow time for the last one to fly back)
    frames(300);
    int stillHere = 0;
    for (const Drone& d : mission.drones()) stillHere += retired.count(d.id) > 0;
    std::printf("  (%zu drones retired, %d not yet removed; moves over fire %d, inside hard zone %d)\n",
                retired.size(), stillHere, crossings, insideHardZone);
    CHECK(stillHere == 0);
    CHECK(crossings == 0);
    CHECK(insideHardZone == 0);

    // --- Grow the fire again: both counts go up on the very next re-plan ---
    const int sweepSmall = mission.activeCount(false), acoSmall = mission.activeCount(true);
    setFire(100);
    replan();
    std::printf("  (fire grows back: sweep %d -> %d, ACO %d -> %d)\n", sweepSmall, mission.activeCount(false),
                acoSmall, mission.activeCount(true));
    CHECK(mission.activeCount(false) > sweepSmall);
    CHECK(mission.activeCount(true) > acoSmall);
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
    CHECK(static_cast<int>(targets.size()) == mission.activeCount(true));

    // ...evenly spaced round it: neighbouring drones are a 1/n share of the tour apart, give or take the
    // longest gap between two waypoints (slots snap to waypoints)
    const std::vector<Vec2>& tour = mission.plan().acoTour;
    const int n = static_cast<int>(tour.size());
    std::vector<float> arc(n, 0.0f);
    float maxGap = 0.0f;
    for (int i = 1; i < n; ++i) {
        const float g = std::sqrt(sqDist(tour[i - 1], tour[i]));
        arc[i] = arc[i - 1] + g;
        maxGap = std::max(maxGap, g);
    }
    const float total = closedRouteLength(tour);
    maxGap = std::max(maxGap, total - arc[n - 1]);
    std::vector<float> pos;
    for (size_t t : targets) pos.push_back(arc[t]);
    const float share = total / pos.size();
    for (size_t i = 0; i < pos.size(); ++i) {
        const float gap = i + 1 < pos.size() ? pos[i + 1] - pos[i] : total - pos[i] + pos[0];
        CHECK(std::fabs(gap - share) <= maxGap + 1e-3f);
    }
}

// After a re-plan, ACO drones already on the tour lock back on at the nearest node (or the next one, if they
// have just passed it) and keep going the same way round: anticlockwise, never back the way they came
TEST(aco_drones_rejoin_at_the_nearest_node_going_forward) {
    const int cx = GRID_WIDTH / 2, cy = GRID_HEIGHT / 2;
    std::vector<std::uint8_t> fire = discFire(cx, cy, 25);
    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
    live.compute(fire);
    Perimeter perim(fire);
    MissionController mission(live, perim.distance(), Vec2(2.5f, GRID_HEIGHT - 2.5f));
    mission.markFireChanged();

    // Anticlockwise on screen = the angle about the fire's centre goes down (y points down the grid)
    auto angle = [&](Vec2 p) { return std::atan2(p.y - (cy + 0.5f), p.x - (cx + 0.5f)); };
    auto turn = [](float from, float to) {   // Signed change of angle, wrapped into (-pi, pi]
        float d = to - from;
        while (d > 3.14159265f) d -= 2.0f * 3.14159265f;
        while (d <= -3.14159265f) d += 2.0f * 3.14159265f;
        return d;
    };

    // Fly long enough for every ACO drone to be out on the tour, then grow the fire so the tour is re-planned
    for (int f = 0; f < 3 * REPLAN_INTERVAL; ++f) mission.update();
    CHECK(loopSignedArea(mission.plan().acoTour) < 0.0);   // Tour handed out anticlockwise
    fire = discFire(cx, cy, 40);
    live.compute(fire);
    perim.update(fire);
    mission.markFireChanged();

    // Run up to and through the re-plan, remembering where each ACO drone was just before it
    std::vector<std::pair<int, Vec2>> before;
    bool replanned = false;
    for (int f = 0; f < REPLAN_INTERVAL + 1 && !replanned; ++f) {
        before.clear();
        for (const Drone& d : mission.drones()) if (d.isACO) before.emplace_back(d.id, d.pos);
        mission.update();
        replanned = mission.takePlanChanged();
    }
    CHECK(replanned);
    const std::vector<Vec2>& tour = mission.plan().acoTour;
    CHECK(loopSignedArea(tour) < 0.0);
    const size_t n = tour.size();

    // Where each drone was sent at the re-plan: the node nearest where it was, or the next one or two (it had
    // passed that node, or reached it during the re-plan frame)
    int checked = 0, wrongNode = 0, backwards = 0;
    std::vector<std::pair<int, float>> followed;   // Drone id, angle at the re-plan
    for (const auto& [id, was] : before) {
        for (const Drone& d : mission.drones()) {
            if (d.id != id || d.state != DroneState::Active) continue;   // Not found: retired on this re-plan
            ++checked;
            const size_t near = nearestWaypoint(tour, was);
            bool ok = false;
            for (size_t k = 0; k <= 2; ++k) ok = ok || d.target == (near + k) % n;
            wrongNode += !ok;
            followed.emplace_back(id, angle(d.pos));
        }
    }

    // Follow them all for a while: each must keep turning anticlockwise round the fire overall
    std::vector<float> total(followed.size(), 0.0f), last(followed.size());
    for (size_t k = 0; k < followed.size(); ++k) last[k] = followed[k].second;
    for (int f = 0; f < 30; ++f) {
        mission.update();
        for (size_t k = 0; k < followed.size(); ++k)
            for (const Drone& d : mission.drones())
                if (d.id == followed[k].first) {
                    total[k] += turn(last[k], angle(d.pos));
                    last[k] = angle(d.pos);
                }
    }
    for (float t : total) backwards += t > 0.0f;
    std::printf("  (%d ACO drones checked: %d sent to another node, %d went backwards)\n", checked, wrongNode, backwards);
    CHECK(checked > 0);
    CHECK(wrongNode == 0);
    CHECK(backwards == 0);
}

TEST(aco_spacing_speed_rule) {
    const float even = 1000.0f;
    CHECK(acoSpacingSpeed(even, even) == 1.0f);                     // Evenly spaced: full speed
    CHECK(acoSpacingSpeed(2.0f * even, even) == 1.0f);              // Big gap ahead: full speed (can't go faster)
    CHECK(acoSpacingSpeed(even * (1.0f - ACO_SPACING_TOLERANCE * 0.5f), even) == 1.0f);   // Within tolerance
    CHECK(std::fabs(acoSpacingSpeed(0.5f * even, even) - 0.5f) < 1e-6f);                  // Half the gap, half speed
    CHECK(acoSpacingSpeed(0.0f, even) == ACO_SPACING_MIN_SPEED);    // Right on the next drone: slowest, never stopped
}

// When a re-plan leaves the ACO drones bunched up, they spread back out to even spacing (without turning back)
TEST(aco_drones_spread_back_out_after_a_replan) {
    const int cx = GRID_WIDTH / 2, cy = GRID_HEIGHT / 2;
    std::vector<std::uint8_t> fire = discFire(cx, cy, 25);
    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
    live.compute(fire);
    Perimeter perim(fire);
    MissionController mission(live, perim.distance(), Vec2(2.5f, GRID_HEIGHT - 2.5f));
    mission.markFireChanged();

    // Worst gap compared with an even share: 1 = perfectly even, 0 = two drones on top of each other
    auto evenness = [&]() {
        const std::vector<float> gaps = mission.acoGaps();
        if (gaps.size() < 2) return 1.0f;
        float total = 0.0f, smallest = std::numeric_limits<float>::max();
        for (float g : gaps) { total += g; smallest = std::min(smallest, g); }
        return smallest / (total / gaps.size());
    };

    // Let the drones settle on the first tour, then bunch them up: move and grow the fire so the tour changes and
    // every drone rejoins where it happens to be.
    // The fire jumps 1.5 x SEARCH_DISTANCE sideways. That puts the whole old tour on the near (west) side of the
    // new one, so every drone rejoins on that side and they end up bunched whatever the starting ACO count in
    // config.h. (A small 30-cell shift only bunched them when there were 5 drones; with ACO_DRONES = 2 the two
    // drones rejoined at 0.85 of an even share and the test failed.) The shift is capped so the new tour still
    // fits on the grid.
    for (int f = 0; f < 6 * REPLAN_INTERVAL; ++f) mission.update();
    const float settled = evenness();
    const int shift = std::min(SEARCH_DISTANCE * 3 / 2, GRID_WIDTH - cx - 45 - SEARCH_DISTANCE - 10);
    fire = discFire(cx + shift, cy, 45);
    live.compute(fire);
    perim.update(fire);
    mission.markFireChanged();
    bool replanned = false;
    for (int f = 0; f <= REPLAN_INTERVAL && !replanned; ++f) {
        mission.update();
        replanned = mission.takePlanChanged();
    }
    CHECK(replanned);
    // Wait until every drone has rejoined the tour, then give the spacing time to even out
    for (int f = 0; f < 2000 && static_cast<int>(mission.acoGaps().size()) < mission.activeCount(true); ++f) mission.update();
    const float rejoined = evenness();
    int slowed = 0, framesToEven = -1;
    for (int f = 0; f < 10 * REPLAN_INTERVAL; ++f) {
        mission.update();
        for (const Drone& d : mission.drones()) slowed += d.speedScale < 1.0f;
        if (framesToEven < 0 && evenness() > 1.0f - 2.0f * ACO_SPACING_TOLERANCE) framesToEven = f + 1;
    }
    const float later = evenness();
    std::printf("  (%d ACO drones; smallest gap / even share: settled %.2f, rejoined after re-plan %.2f, even again "
                "(> %.2f) after %.0f simulated min, %.2f after %.0f min)\n", mission.activeCount(true), settled,
                rejoined, 1.0f - 2.0f * ACO_SPACING_TOLERANCE, framesToEven * SIM_MINUTES_PER_TICK, later,
                10 * REPLAN_INTERVAL * SIM_MINUTES_PER_TICK);
    CHECK(settled > 1.0f - 2.0f * ACO_SPACING_TOLERANCE);    // Even to begin with
    CHECK(rejoined < 0.8f);                                  // The re-plan really did bunch them up...
    CHECK(slowed > 0);                                       // ...the crowded ones slowed down...
    CHECK(later > 1.0f - 2.0f * ACO_SPACING_TOLERANCE);      // ...and the spacing is even again
    // Once even, everybody is back to full speed
    for (const Drone& d : mission.drones()) CHECK(d.speedScale == 1.0f);
}

// The ACO count follows the one-drone loop time: each point is revisited within ACO_REVISIT_S
TEST(aco_count_follows_the_revisit_time) {
    for (int r : { 1, 100, 200 }) {
        SweepSetup s(discFire(CX, CY, r));
        const MissionPlan& p = s.plan;
        std::printf("  (fire %d m across: tour %.1f km, %.1f min for one drone, %d ACO drones)\n",
                    2 * r * static_cast<int>(SIM_CELL_SIZE_M), closedRouteLength(p.acoTour) * SIM_CELL_SIZE_M / 1000.0f,
                    p.acoLoopS / 60.0f, p.acoNeeded);
        CHECK(std::fabs(p.acoLoopS - flightSeconds(closedRouteLength(p.acoTour))) < 1e-3f);
        CHECK(p.acoNeeded >= MIN_ACO_DRONES && p.acoNeeded <= MAX_ACO_DRONES);
        CHECK(p.acoLoopS / p.acoNeeded <= ACO_REVISIT_S);                    // Revisit within the limit...
        CHECK(p.acoNeeded == MIN_ACO_DRONES || p.acoLoopS / (p.acoNeeded - 1) > ACO_REVISIT_S);   // ...with no spare
        CHECK(!p.acoOverLimit);
    }
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
