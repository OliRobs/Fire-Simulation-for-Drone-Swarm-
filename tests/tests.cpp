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
#include <queue>
#include <random>
#include <set>
#include <string>
#include <unordered_map>
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
#include "wind_overlay.h"
#include "wind_solver.h"

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

// The live distance field is now updated incrementally (only cells that caught fire since the last frame are
// added) and clamped: exact out to maxDistance + 1, "very far" beyond. At every step it must equal the exact
// range-limited transform with everything beyond that clamp radius read as very far. Covers a real fire growing
// frame by frame (plus a second fire started part-way, which grows the window), and the fallbacks to a full
// transform: the first call, a source cell switching off, and a burst too big for per-cell updates.
TEST(live_edt_incremental_matches_clamped_full_transform) {
    const int range = LIVE_EDT_RANGE;
    const double cap = (range + 1.0) * (range + 1.0);
    DistanceTransform inc(GRID_WIDTH, GRID_HEIGHT, range, DistanceTransform::Mode::IncrementalClamped);
    DistanceTransform exact(GRID_WIDTH, GRID_HEIGHT, range);
    int steps = 0, wrong = 0, windowWrong = 0;
    std::uint64_t lastVersion = inc.version();
    // Compares the two fields everywhere either has a finite value (outside both windows both read very far)
    auto compare = [&](const std::vector<std::uint8_t>& mask, const GridRect& bounds) {
        inc.compute(mask, bounds);
        exact.compute(mask, bounds);
        const GridRect w = inc.window().unite(exact.window());
        for (int y = w.y0; y < w.y1; ++y)
            for (int x = w.x0; x < w.x1; ++x) {
                const double e = exact.squaredDistance(at(x, y)), v = inc.squaredDistance(at(x, y));
                wrong += e <= cap ? v != e : !(v > cap);
            }
        const GridRect a = inc.window(), b = exact.window();
        windowWrong += a.x0 != b.x0 || a.y0 != b.y0 || a.x1 != b.x1 || a.y1 != b.y1;
        wrong += inc.version() == lastVersion;   // Every compute() is a new version
        lastVersion = inc.version();
        ++steps;
    };

    // --- A real fire, frame by frame, with a second fire part-way ---
    std::vector<Cell> grid(N);
    FireSimulation sim(grid);
    sim.ignite(at(CX, CY));
    for (int t = 0; t < 500; ++t) {
        sim.step();
        if (t == 200) sim.ignite(at(CX + 150, CY - 100));
        if (sim.takeFireGrew()) compare(sim.fireMask(), sim.fireBounds());
    }
    const int growthSteps = steps;

    // --- Fallbacks ---
    std::vector<std::uint8_t> mask = sim.fireMask();
    GridRect box = sim.fireBounds();
    // A source switches off (as if the fire were reset there): distances grow, so a full transform is needed
    for (int x = box.x0; x < box.x1; ++x) mask[at(x, (box.y0 + box.y1) / 2)] = 0;
    compare(mask, box);
    // A burst: a whole disc catches fire in one step (as a crown-fire run might)
    const auto burst = discFire(CX - 100, CY + 80, 90);
    for (int i = 0; i < N; ++i) if (burst[i]) { mask[i] = 1; box.include(i % GRID_WIDTH, i / GRID_WIDTH); }
    compare(mask, box);
    // Back to small additions after the fallbacks
    for (int k = 0; k < 5; ++k) { mask[at(box.x0 + k, box.y0)] = 1; compare(mask, box); }

    std::printf("  (%d growth steps and %d fallback / follow-up steps compared)\n", growthSteps, steps - growthSteps);
    CHECK(growthSteps > 100);
    CHECK(wrong == 0);
    CHECK(windowWrong == 0);
}

// The incremental live field lowers distances one source at a time with DistanceTransform::addSource, whose
// inner loop was rewritten so the compiler can vectorise it, with the disc's row half-widths precomputed. This
// checks that path against a brute-force answer (minimum squared distance to every source, "very far" beyond
// the clamp radius) where the real-fire test above cannot reach: discs clipped by every grid edge and corner,
// sources on the outermost rows and columns, grids narrower than the disc (down to one column), and several clamp radii (1 cell, a
// small one, and the live field's own LIVE_EDT_RANGE). Sources are added one per compute() so the field takes
// the per-cell route rather than the full-transform fallback.
TEST(live_edt_add_source_matches_brute_force_at_edges) {
    struct Case { int w, h, range; };
    long long compared = 0, wrong = 0;
    std::mt19937 rng(20261127);
    for (const Case c : { Case{ 200, 150, LIVE_EDT_RANGE }, Case{ 45, 300, LIVE_EDT_RANGE }, Case{ 97, 61, 7 },
                          Case{ 64, 64, 1 }, Case{ 1, 120, 1 } }) {
        const double cap = (c.range + 1.0) * (c.range + 1.0);
        DistanceTransform inc(c.w, c.h, c.range, DistanceTransform::Mode::IncrementalClamped);
        std::vector<std::uint8_t> mask(static_cast<size_t>(c.w) * c.h, 0);
        GridRect box;
        std::vector<std::pair<int, int>> sources;

        // Corners, edge midpoints and the cells next to them, then random cells (some on the border rows)
        std::vector<std::pair<int, int>> order = {
            { 0, 0 }, { c.w - 1, c.h - 1 }, { c.w - 1, 0 }, { 0, c.h - 1 }, { c.w / 2, 0 }, { c.w / 2, c.h - 1 },
            { 0, c.h / 2 }, { c.w - 1, c.h / 2 }, { std::min(1, c.w - 1), 1 }, { std::max(0, c.w - 2), c.h - 2 } };
        std::uniform_int_distribution<int> rx(0, c.w - 1), ry(0, c.h - 1);
        for (int k = 0; k < 30; ++k) order.push_back({ rx(rng), k % 5 == 0 ? 0 : ry(rng) });

        for (const auto& s : order) {
            const size_t idx = static_cast<size_t>(s.second) * c.w + s.first;
            if (mask[idx]) continue;   // Already a source (random picks can repeat)
            mask[idx] = 1;
            box.include(s.first, s.second);
            sources.push_back(s);
            inc.compute(mask, box);

            // Brute force: the exact squared distance to the nearest source, everywhere on the grid
            for (int y = 0; y < c.h; ++y)
                for (int x = 0; x < c.w; ++x) {
                    double e = std::numeric_limits<double>::infinity();
                    for (const auto& q : sources) {
                        const double dx = x - q.first, dy = y - q.second;
                        e = std::min(e, dx * dx + dy * dy);
                    }
                    const double v = inc.squaredDistance(y * c.w + x);
                    wrong += e <= cap ? v != e : !(v > cap);   // Exact within the clamp, "very far" beyond it
                    ++compared;
                }
        }
    }
    std::printf("  (%lld cell values compared)\n", compared);
    CHECK(compared > 1000000);
    CHECK(wrong == 0);
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

// The perimeter update now tests only the cells that caught fire since the last call (live cells already
// tested can't trigger again: see Perimeter::update). It must move exactly the same sections as the original,
// which tested every live cell of the fire's area on every call. This keeps the original word for word and
// compares the two on a real growing fire (with a second fire lit part-way) and on a mask that is not
// growing-only: discs that jump about, so cells go out and catch fire again, some of them far out.
namespace {
    class OriginalPerimeter {
    public:
        explicit OriginalPerimeter(const std::vector<std::uint8_t>& initialFire)
            : committed_(initialFire), trigger_(initialFire.size(), 0),
            perimEdt_(GRID_WIDTH, GRID_HEIGHT, SEARCH_DISTANCE + 1),
            triggerEdt_(GRID_WIDTH, GRID_HEIGHT, PERIMETER_SECTION_RADIUS + 1) {
            perimEdt_.compute(committed_);
            int x0, y0, x1, y1;
            if (maskBounds(committed_, GRID_WIDTH, GRID_HEIGHT, x0, y0, x1, y1)) fireSeen_ = GridRect{ x0, y0, x1 + 1, y1 + 1 };
        }
        // The update as it was: every live cell of the area tested against the trigger distance
        bool update(const std::vector<std::uint8_t>& liveFire, const GridRect& fireBounds) {
            constexpr double triggerSq = static_cast<double>(PERIMETER_DISTANCE - PERIMETER_UPDATE_DISTANCE)
                * (PERIMETER_DISTANCE - PERIMETER_UPDATE_DISTANCE);
            constexpr double sectionSq = static_cast<double>(PERIMETER_SECTION_RADIUS) * PERIMETER_SECTION_RADIUS;
            fireSeen_ = fireSeen_.unite(fireBounds);
            const GridRect area = fireSeen_;
            bool any = false;
            for (int y = area.y0; y < area.y1; ++y) {
                for (int idx = y * GRID_WIDTH + area.x0, end = y * GRID_WIDTH + area.x1; idx < end; ++idx) {
                    const bool trig = liveFire[idx] && perimEdt_.squaredDistance(idx) >= triggerSq;
                    trigger_[idx] = trig;
                    any = any || trig;
                }
            }
            if (!any) return false;
            triggerEdt_.compute(trigger_, area);
            for (int y = area.y0; y < area.y1; ++y) {
                for (int idx = y * GRID_WIDTH + area.x0, end = y * GRID_WIDTH + area.x1; idx < end; ++idx) {
                    if (liveFire[idx] && !committed_[idx] && triggerEdt_.squaredDistance(idx) <= sectionSq) committed_[idx] = 1;
                }
            }
            perimEdt_.compute(committed_, area);
            return true;
        }
        const DistanceTransform& distance() const { return perimEdt_; }
    private:
        std::vector<std::uint8_t> committed_, trigger_;
        DistanceTransform perimEdt_, triggerEdt_;
        GridRect fireSeen_;
    };

    // Bounding box of a mask's non-zero cells (empty if none), as FireSimulation::fireBounds() gives it
    GridRect boundsOf(const std::vector<std::uint8_t>& m) {
        int x0, y0, x1, y1;
        if (!maskBounds(m, GRID_WIDTH, GRID_HEIGHT, x0, y0, x1, y1)) return GridRect{};
        return GridRect{ x0, y0, x1 + 1, y1 + 1 };
    }

    // True if two committed-perimeter fields are identical, bit for bit, at every cell of the grid
    bool sameField(const DistanceTransform& a, const DistanceTransform& b) {
        for (int i = 0; i < N; ++i) {
            const double da = a.squaredDistance(i), db = b.squaredDistance(i);
            if (!(da == db)) return false;   // (INF == INF holds; written this way so a NaN would fail too)
        }
        return true;
    }
}

TEST(perimeter_update_tests_only_new_cells_and_matches_original) {
    int calls = 0, moves = 0, mismatches = 0;

    // --- A real fire: one ignition, then a second, separate one part-way through ---
    {
        std::vector<Cell> grid(N);
        FireSimulation sim(grid);
        sim.ignite(at(CX, CY));
        for (int t = 0; t < 100; ++t) sim.step();
        Perimeter fast(sim.fireMask());
        OriginalPerimeter original(sim.fireMask());
        for (int t = 0; t < 2400; ++t) {
            if (t == 900) sim.ignite(at(CX + 150, CY + 80));   // A spot fire well outside the perimeter
            sim.step();
            if (!sim.takeFireGrew()) continue;
            const bool a = fast.update(sim.fireMask(), sim.fireBounds());
            const bool b = original.update(sim.fireMask(), sim.fireBounds());
            ++calls;
            mismatches += a != b;
            if (a) {
                ++moves;
                mismatches += !sameField(fast.distance(), original.distance());
            }
        }
        mismatches += !sameField(fast.distance(), original.distance());
    }

    // --- A mask that is not growing-only: discs that move about and change size, so cells go out and catch
    //     fire again, including cells far beyond the trigger distance that were live before ---
    {
        std::mt19937 rng(7);
        std::uniform_int_distribution<int> dx(-60, 60), rad(5, 45);
        std::vector<std::uint8_t> mask = discFire(CX, CY, 20);
        Perimeter fast(mask);
        OriginalPerimeter original(mask);
        for (int step = 0; step < 60; ++step) {
            mask = discFire(CX + dx(rng), CY + dx(rng), rad(rng));
            if (step % 3 == 0) {   // Sometimes add a second, distant disc as well
                const auto other = discFire(CX + 4 * dx(rng), CY + 3 * dx(rng), rad(rng) / 2);
                for (int i = 0; i < N; ++i) mask[i] = mask[i] | other[i];
            }
            const GridRect b = boundsOf(mask);
            const bool ra = fast.update(mask, b);
            const bool rb = original.update(mask, b);
            ++calls;
            mismatches += ra != rb;
            if (ra) ++moves;
            if (ra || rb) mismatches += !sameField(fast.distance(), original.distance());
        }
        mismatches += !sameField(fast.distance(), original.distance());
    }

    std::printf("  (%d updates compared, %d perimeter moves, %d differences)\n", calls, moves, mismatches);
    CHECK(moves >= 5);         // The perimeter really moved, many times
    CHECK(mismatches == 0);
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

// ----------------------------------------------------------------------------
// Original (pre-2026-11-20) SafePlanner, kept word for word apart from the class wrapper, as the reference for
// the bitmap / bucketed-open-list / one-record-per-cell A*. Paths must match it point for point.
// ----------------------------------------------------------------------------
namespace ref {
    class OriginalPlanner {
    public:
        explicit OriginalPlanner(const DistanceTransform& liveEdt)
            : edt_(liveEdt), g_(GRID_WIDTH * GRID_HEIGHT), came_(GRID_WIDTH * GRID_HEIGHT),
              openStamp_(GRID_WIDTH * GRID_HEIGHT, 0), closedStamp_(GRID_WIDTH * GRID_HEIGHT, 0) {}
        bool cellOk(int x, int y) const;
        bool lineClear(Vec2 a, Vec2 b) const;
        std::vector<Vec2> leg(Vec2 a, Vec2 b);
    private:
        int snap(Vec2 p) const;
        float heuristic(int a, int b) const;
        const DistanceTransform& edt_;
        std::vector<float> g_;
        std::vector<int> came_;
        std::vector<std::uint32_t> openStamp_, closedStamp_;
        std::uint32_t search_ = 0;
    };

bool OriginalPlanner::cellOk(int x, int y) const {
    return x >= 0 && x < GRID_WIDTH && y >= 0 && y < GRID_HEIGHT
        && edt_.squaredDistance(y * GRID_WIDTH + x) >= HARD_NO_FLY_SQ;
}

bool OriginalPlanner::lineClear(Vec2 a, Vec2 b) const {
    const float len = std::sqrt(sqDist(a, b));
    const int steps = std::max(1, static_cast<int>(len / 0.25f) + 1);

    for (int i = 0; i <= steps; ++i) {
        const float t = static_cast<float>(i) / steps;
        const Vec2 p = a + (b - a) * t;
        if (!cellOk(static_cast<int>(std::floor(p.x)), static_cast<int>(std::floor(p.y)))) {
            return false;
        }
    }
    return true;
}

std::vector<Vec2> OriginalPlanner::leg(Vec2 a, Vec2 b) {
    // Fast path: direct line of sight is clear
    if (lineClear(a, b)) return { b };

    // Snap start and destination to the nearest flyable cells. If either has no flyable cell
    // nearby there is nothing safe to fly to, so report failure rather than a straight line.
    const int s = snap(a), t = snap(b);
    if (s < 0 || t < 0) return {};

    // Does b itself lie in flyable space? If not, the path will finish at the snapped cell instead.
    const bool goalFlyable = cellOk(static_cast<int>(std::floor(b.x)), static_cast<int>(std::floor(b.y)));

    // --- Phase 1: A* Grid Pathfinding ---
    // New search: bumping the stamp invalidates every g_/came_ entry from previous searches at once.
    if (++search_ == 0) {   // Stamp wrapped round (after ~4 billion searches): reset properly
        std::fill(openStamp_.begin(), openStamp_.end(), 0);
        std::fill(closedStamp_.begin(), closedStamp_.end(), 0);
        search_ = 1;
    }
    auto gOf = [&](int v) {   // Cost so far, or infinity if v has not been reached in this search
        return openStamp_[v] == search_ ? g_[v] : std::numeric_limits<float>::max();
    };

    using QItem = std::pair<float, int>; // Priority queue item: {f_score, cell_index}
    std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> open;

    g_[s] = 0.0f;
    came_[s] = -1;
    openStamp_[s] = search_;
    open.push({ heuristic(s, t), s });

    bool reached = false;
    while (!open.empty()) {
        const int u = open.top().second;
        open.pop();

        if (closedStamp_[u] == search_) continue; // Stale queue entry: already expanded
        closedStamp_[u] = search_;

        if (u == t) { reached = true; break; }    // Goal reached

        const int ux = u % GRID_WIDTH, uy = u / GRID_WIDTH;

        // Expand all 8 neighbours
        for (int d = 0; d < NUM_DIRS; ++d) {
            const int nx = ux + DX[d], ny = uy + DY[d];
            if (!cellOk(nx, ny)) continue;

            // No corner cutting: a diagonal step needs both orthogonal cells it squeezes past to be flyable
            const bool diagonal = DX[d] != 0 && DY[d] != 0;
            if (diagonal && (!cellOk(ux + DX[d], uy) || !cellOk(ux, uy + DY[d]))) continue;

            const int v = ny * GRID_WIDTH + nx;
            const float ng = g_[u] + (diagonal ? 1.41421f : 1.0f);

            if (ng < gOf(v)) {
                g_[v] = ng;
                came_[v] = u;
                openStamp_[v] = search_;
                open.push({ ng + heuristic(v, t), v });
            }
        }
    }

    // The goal is cut off from the start by fire: there is no safe leg.
    if (!reached) return {};

    // --- Phase 2: Reconstruct the cell path (start point, cell centres..., end point) ---
    std::vector<Vec2> full;
    for (int v = t; v != -1; v = came_[v]) {
        full.push_back(cellCentre(v));
    }
    full.push_back(a);
    std::reverse(full.begin(), full.end());
    if (goalFlyable) full.push_back(b);   // Otherwise finish at the centre of the snapped cell

    // --- Phase 3: String pulling ---
    // From each point, jump to the farthest later point that is still in clear line of sight.
    std::vector<Vec2> out;
    size_t i = 0;
    while (i + 1 < full.size()) {
        size_t j = full.size() - 1;
        while (j > i + 1 && !lineClear(full[i], full[j])) {
            --j;
        }
        out.push_back(full[j]);
        i = j;
    }
    return out;
}

int OriginalPlanner::snap(Vec2 p) const {
    const int cx = static_cast<int>(std::floor(p.x)), cy = static_cast<int>(std::floor(p.y));

    for (int r = 0; r <= 40; ++r) {
        for (int dy = -r; dy <= r; ++dy) {
            for (int dx = -r; dx <= r; ++dx) {
                if (std::max(std::abs(dx), std::abs(dy)) == r && cellOk(cx + dx, cy + dy)) {
                    return (cy + dy) * GRID_WIDTH + (cx + dx);
                }
            }
        }
    }
    return -1; // No valid cell found within search radius
}

float OriginalPlanner::heuristic(int a, int b) const {
    const float dx = static_cast<float>(std::abs(a % GRID_WIDTH - b % GRID_WIDTH));
    const float dy = static_cast<float>(std::abs(a / GRID_WIDTH - b / GRID_WIDTH));
    return dx + dy - 0.58579f * std::min(dx, dy);
}
}   // namespace ref

// A* now reads a flyable-cell bitmap, keeps its open list in f-cost buckets and keeps each cell's state in one
// record. It must return exactly the legs of the original planner (ref::OriginalPlanner), point for point: on
// several fire shapes, with the range-limited live transform the app uses and the whole-grid one, including legs
// whose goal is inside the no-fly zone or cut off, and after the fire grows (the bitmap must follow the field).
TEST(astar_matches_original_planner) {
    std::vector<std::uint8_t> twoFires = discFire(CX - 120, CY, 60);
    const auto other = discFire(CX + 90, CY + 40, 45);
    for (int i = 0; i < N; ++i) twoFires[i] |= other[i];
    std::vector<std::uint8_t> edgeFire = discFire(15, CY, 70);
    std::mt19937 rng(5);
    std::uniform_real_distribution<float> px(1.0f, GRID_WIDTH - 1.0f), py(1.0f, GRID_HEIGHT - 1.0f);
    int legs = 0, detours = 0, empty = 0, different = 0;
    for (int range : { LIVE_EDT_RANGE, 0 }) {
        DistanceTransform live(GRID_WIDTH, GRID_HEIGHT, range);
        SafePlanner fast(live);
        ref::OriginalPlanner original(live);
        // The same planners carry on through every fire, so the bitmap has to be refreshed as the field changes
        for (const auto* fire : { &twoFires, &edgeFire }) {
            for (int grow = 0; grow < 2; ++grow) {
                std::vector<std::uint8_t> f = *fire;
                if (grow) { const auto more = discFire(CX, CY - 150, 50); for (int i = 0; i < N; ++i) f[i] |= more[i]; }
                live.compute(f);
                for (int k = 0; k < 60; ++k) {
                    // Every 10th goal is the centre of a fire (inside the no-fly zone)
                    const Vec2 a(px(rng), py(rng));
                    const Vec2 b = k % 10 == 0 ? Vec2(CX - 120.5f, CY + 0.5f) : Vec2(px(rng), py(rng));
                    const std::vector<Vec2> p = fast.leg(a, b), q = original.leg(a, b);
                    bool same = p.size() == q.size();
                    for (size_t i = 0; same && i < p.size(); ++i) same = p[i].x == q[i].x && p[i].y == q[i].y;
                    different += !same;
                    detours += q.size() > 1;
                    empty += q.empty();
                    ++legs;
                }
            }
        }
    }
    std::printf("  (%d legs, %d of them A* detours, %d with no safe path; %d differ)\n", legs, detours, empty, different);
    CHECK(different == 0);
    CHECK(detours > 50);   // The legs really exercise A*
}

// lineClear skips samples that the distance field proves are clear (see SafePlanner::lineClear). It must give
// exactly the answer of the original loop that tests every sample (ref::OriginalPlanner::lineClear), for every
// kind of field a planner can be given: the app's incremental clamped live field (also after the fire grows, so
// the incremental path runs), a range-limited exact one, the whole-grid one, and one with no fire at all. The
// segments include random ones across the map, ones that graze the no-fly radius round a fire (where a wrong skip
// would show), and ones touching the grid edges (where the skip is turned off).
TEST(lineclear_skip_matches_every_sample) {
    const int R = 50;                                  // Fire disc radius (cells)
    std::vector<std::uint8_t> fire = discFire(CX, CY, R);
    std::vector<std::uint8_t> grown = fire;
    const auto more = discFire(CX + 70, CY - 30, 35);
    for (int i = 0; i < N; ++i) grown[i] |= more[i];
    const std::vector<std::uint8_t> none(N, 0);

    std::mt19937 rng(11);
    std::uniform_real_distribution<float> px(0.0f, GRID_WIDTH - 0.001f), py(0.0f, GRID_HEIGHT - 0.001f);
    std::uniform_real_distribution<float> ang(0.0f, 6.2831853f), near(-3.0f, 3.0f), u(0.0f, 1.0f);
    const float hard = static_cast<float>(std::sqrt(HARD_NO_FLY_SQ));
    int segments = 0, clear = 0, blocked = 0, differ = 0;

    // Compares the two planners on one segment and keeps count
    auto compare = [&](const SafePlanner& fast, const ref::OriginalPlanner& orig, Vec2 a, Vec2 b) {
        const bool f = fast.lineClear(a, b), o = orig.lineClear(a, b);
        differ += f != o;
        (o ? clear : blocked) += 1;
        ++segments;
    };

    struct FieldCase { int range; DistanceTransform::Mode mode; };
    for (const FieldCase fc : { FieldCase{ LIVE_EDT_RANGE, DistanceTransform::Mode::IncrementalClamped },
                                FieldCase{ LIVE_EDT_RANGE, DistanceTransform::Mode::Exact },
                                FieldCase{ 0, DistanceTransform::Mode::Exact } }) {
        DistanceTransform live(GRID_WIDTH, GRID_HEIGHT, fc.range, fc.mode);
        const SafePlanner fast(live);
        const ref::OriginalPlanner orig(live);
        const std::vector<std::uint8_t>* fields[] = { &none, &fire, &grown };
        for (const std::vector<std::uint8_t>* f : fields) {   // Same field object throughout: incremental after the first
            live.compute(*f);
            for (int k = 0; k < 400; ++k) {
                // Random segment anywhere on the map
                compare(fast, orig, Vec2(px(rng), py(rng)), Vec2(px(rng), py(rng)));
                // Grazing: a chord of a circle whose radius is within a few cells of the fire disc plus the no-fly
                // radius, so its closest approach to the fire is right at the limit
                const float rad = R + hard + near(rng), a0 = ang(rng), a1 = a0 + 0.3f + u(rng) * 2.0f;
                const float c = rad / std::cos((a1 - a0) / 2.0f);   // Endpoints so the chord's midpoint is at rad
                compare(fast, orig, Vec2(CX + 0.5f + c * std::cos(a0), CY + 0.5f + c * std::sin(a0)),
                                    Vec2(CX + 0.5f + c * std::cos(a1), CY + 0.5f + c * std::sin(a1)));
                // Along or close to a grid edge
                compare(fast, orig, Vec2(px(rng), u(rng) * 2.0f), Vec2(px(rng), GRID_HEIGHT - u(rng) * 2.0f));
                compare(fast, orig, Vec2(u(rng), py(rng)), Vec2(GRID_WIDTH - u(rng), py(rng)));
            }
        }
    }
    std::printf("  (%d segments: %d clear, %d blocked; %d differ)\n", segments, clear, blocked, differ);
    CHECK(differ == 0);
    CHECK(blocked > 1000);   // Plenty of both answers, so the comparison means something
    CHECK(clear > 1000);
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

// --- Rothermel set-up: cached fuel-bed terms and the per-cell memos give bit-identical rates ---
namespace {
    // The wind and slope factors and the effective wind speed exactly as they were written before the fuel-bed
    // terms were cached (kept word for word, so the new code can be checked against them bit for bit)
    double origWindFactor(const SurfaceFuelBed& b, double u) {
        if (u <= 0.0) return 0.0;
        return b.windC * std::pow(u, b.windB) * std::pow(b.relativePacking, -b.windE);
    }
    double origSlopeFactor(const SurfaceFuelBed& b, double tanSlope) {
        if (b.packingRatio <= 0.0) return 0.0;
        return 5.275 * std::pow(b.packingRatio, -0.3) * tanSlope * tanSlope;
    }
    double origEffectiveWindSpeed(const SurfaceFuelBed& b, double phiE) {
        if (phiE <= 0.0 || b.windC <= 0.0 || b.windB <= 0.0) return 0.0;
        return std::pow(phiE * std::pow(b.relativePacking, b.windE) / b.windC, 1.0 / b.windB);
    }

    // The original computeSpread, word for word, on the original factors above
    SpreadEllipse origComputeSpread(const SurfaceFuelBed& bed, double midflame, double windX, double windY,
                                    double tanSlope, double upX, double upY, bool applyWindLimit) {
        SpreadEllipse e;
        if (!bed.burnable()) return e;
        const double phiW = origWindFactor(bed, midflame);
        const double phiS = origSlopeFactor(bed, tanSlope);
        const double vx = phiS * upX + phiW * windX;
        const double vy = phiS * upY + phiW * windY;
        double phiE = std::sqrt(vx * vx + vy * vy);
        if (phiE > 1e-12) { e.dirX = vx / phiE; e.dirY = vy / phiE; }
        double uE = origEffectiveWindSpeed(bed, phiE);
        if (applyWindLimit) {
            const double uMax = 0.9 * bed.reactionIntensity;
            if (uE > uMax) {
                uE = uMax;
                phiE = origWindFactor(bed, uE);
            }
        }
        e.headRos = bed.r0 * (1.0 + phiE);
        e.effectiveWindMph = uE / FT_PER_MIN_PER_MPH;
        const double z = 1.0 + 0.25 * e.effectiveWindMph;
        e.eccentricity = std::sqrt(z * z - 1.0) / z;
        return e;
    }

    // Moisture scenarios used by both tests: the D2L2 default, very dry, fully cured, and wet enough that some
    // models stop burning (so non-burnable beds, with all-zero fields, are exercised too)
    std::vector<FuelMoisture> testMoistures() {
        FuelMoisture dry, cured, wet;
        dry.dead1h = 0.03; dry.dead10h = 0.04; dry.dead100h = 0.05; dry.liveHerb = 0.30; dry.liveWoody = 0.60;
        cured.liveHerb = 0.30;
        wet.dead1h = 0.20; wet.dead10h = 0.22; wet.dead100h = 0.24; wet.liveHerb = 1.50; wet.liveWoody = 1.50;
        return { FuelMoisture(), dry, cured, wet };
    }
}

// The wind factor, slope factor, effective wind speed and whole spread ellipse, now using the fuel-bed terms
// cached by computeFuelBed, must equal the original inline-pow formulas bit for bit: every fuel model (plus
// fuel model 10), four moisture scenarios, and a spread of wind speeds, slopes and directions, with and
// without the wind limit. computeSpreadWithWindFactor with a precomputed phi_w must equal computeSpread.
TEST(rothermel_cached_bed_terms_match_original_formulas) {
    int count = 0;
    const FuelModel* models = allFuelModels(count);
    std::vector<const FuelModel*> list;
    for (int i = 0; i < count; ++i) list.push_back(&models[i]);
    list.push_back(&fuelModel10());

    std::mt19937 rng(11);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    long compared = 0, differ = 0;
    int burnableBeds = 0, deadBeds = 0;
    for (const FuelMoisture& mo : testMoistures()) {
        for (const FuelModel* m : list) {
            const SurfaceFuelBed bed = computeFuelBed(*m, mo);
            (bed.burnable() ? burnableBeds : deadBeds)++;
            for (int k = 0; k < 200; ++k) {
                // Wind 0..60 mi/h midflame (k = 0: exactly 0), slope 0..45 degrees, random wind and aspect
                const double u = k == 0 ? 0.0 : 60.0 * FT_PER_MIN_PER_MPH * unit(rng);
                const double tanSlope = k % 7 == 0 ? 0.0 : std::tan(0.785 * unit(rng));
                const double phiE = 50.0 * unit(rng);
                const double wa = 6.283185307 * unit(rng), sa = 6.283185307 * unit(rng);
                const bool limit = k % 2 == 1;
                ++compared;
                if (bed.windFactor(u) != origWindFactor(bed, u)) ++differ;
                if (bed.slopeFactor(tanSlope) != origSlopeFactor(bed, tanSlope)) ++differ;
                if (bed.effectiveWindSpeed(phiE) != origEffectiveWindSpeed(bed, phiE)) ++differ;
                const SpreadEllipse a = computeSpread(bed, u, std::cos(wa), std::sin(wa), tanSlope,
                                                      std::cos(sa), std::sin(sa), limit);
                const SpreadEllipse b = origComputeSpread(bed, u, std::cos(wa), std::sin(wa), tanSlope,
                                                          std::cos(sa), std::sin(sa), limit);
                const SpreadEllipse c = computeSpreadWithWindFactor(bed, bed.windFactor(u), std::cos(wa),
                                                                    std::sin(wa), tanSlope, std::cos(sa),
                                                                    std::sin(sa), limit);
                for (const SpreadEllipse* s : { &a, &c }) {
                    if (s->headRos != b.headRos || s->dirX != b.dirX || s->dirY != b.dirY
                        || s->eccentricity != b.eccentricity || s->effectiveWindMph != b.effectiveWindMph) ++differ;
                }
            }
        }
    }
    std::printf("  %ld cases (%d burnable beds, %d that won't burn), %ld differ\n", compared, burnableBeds,
                deadBeds, differ);
    CHECK(differ == 0);
    CHECK(burnableBeds > 0 && deadBeds > 0);   // Both kinds of bed were exercised
}

// FireSimulation's set-up (fuel-bed terms cached, wind factor memo per fuel model, unsheltered WAF per fuel
// model, crown wind factor computed once) must give exactly the per-cell rates of the original set-up loop,
// kept word for word below on the original formulas. The grid mixes every fuel model in patches, random
// slopes and aspects, forest patches with crown fire (active and passive), awkward canopy edge cases (cover
// without height, height without cover, sparse cover on the 5 % crown fill boundary, canopy with no bulk
// density), roads, non-burnable
// codes and code 0. Several weathers: the default, a fixed WAF, the wind limit, no crown fire, calm, and a
// wet scenario. Cells are 10 x 12 m so the step directions are not the square-grid ones.
TEST(fire_setup_matches_original_rates) {
    int count = 0;
    const FuelModel* models = allFuelModels(count);
    std::mt19937 rng(23);
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    std::vector<Cell> grid(N);
    for (int y = 0; y < GRID_HEIGHT; ++y) {
        for (int x = 0; x < GRID_WIDTH; ++x) {
            Cell& c = grid[at(x, y)];
            const int patch = (x / 30) * 7 + (y / 30) * 13;
            c.fuelModel = static_cast<std::int16_t>(models[patch % count].code);
            if (patch % 53 == 0) c.fuelModel = 91;   // Non-burnable (urban)
            if (patch % 61 == 0) c.fuelModel = 0;    // No data
            if (x % 97 == 0) c.roadClass = static_cast<std::uint8_t>(20 + y % 4);   // A road: fuel break
            c.slopeDeg = (x / 30 + y / 30) % 5 == 0 ? 0.0f : 45.0f * unit(rng);   // Some flat patches
            const float a = 6.2831853f * unit(rng);
            c.upslopeX = c.slopeDeg > 0.0f ? std::cos(a) : 0.0f;
            c.upslopeY = c.slopeDeg > 0.0f ? std::sin(a) : 0.0f;
            switch (((x / 40) + 3 * (y / 40)) % 6) {
            case 0: case 1:   // Forest
                c.canopyCoverPct = 10.0f + 80.0f * unit(rng);
                c.canopyHeightM = 5.0f + 30.0f * unit(rng);
                c.canopyBaseHeightM = 0.3f + 6.0f * unit(rng);
                c.canopyBulkDensity = 0.02f + 0.35f * unit(rng);
                break;
            case 2:           // Odd canopy layers: cover without height, height without cover, sparse cover,
                              // and a canopy with no bulk density (shelters the wind, can't crown)
                if (x % 4 == 0) c.canopyCoverPct = 50.0f;
                else if (x % 4 == 1) { c.canopyHeightM = 20.0f; c.canopyBulkDensity = 0.1f; }
                else if (x % 4 == 2) { c.canopyCoverPct = 15.0f * unit(rng); c.canopyHeightM = 15.0f;
                                       c.canopyBaseHeightM = 15.0f * unit(rng); c.canopyBulkDensity = 0.1f; }
                else { c.canopyCoverPct = 60.0f; c.canopyHeightM = 20.0f; c.canopyBaseHeightM = 2.0f; }
                break;
            default: break;   // Open ground: no canopy
            }
        }
    }

    std::vector<FireWeather> weathers;
    weathers.push_back(FireWeather());
    { FireWeather w; w.windAdjustmentFactor = 0.4f; weathers.push_back(w); }
    { FireWeather w; w.windLimit = !w.windLimit; w.wind20ftKmh = 60.0f; weathers.push_back(w); }
    { FireWeather w; w.crownFire = false; w.windFromDeg = 135.0f; weathers.push_back(w); }
    { FireWeather w; w.wind20ftKmh = 0.0f; weathers.push_back(w); }
    { FireWeather w; w.moisture = testMoistures()[3]; w.wind20ftKmh = 45.0f; weathers.push_back(w); }

    const float csx = 10.0f, csy = 12.0f;
    long cellsCompared = 0, differ = 0, active = 0, passive = 0, sheltered = 0;
    for (const FireWeather& weather : weathers) {
        FireSimulation sim(grid, weather, csx, csy);

        // --- The original set-up loop, word for word except for writing to local arrays ---
        double stepUx[NUM_DIRS], stepUy[NUM_DIRS];
        for (int d = 0; d < NUM_DIRS; ++d) {
            const double sx = DX[d] * csx, sy = DY[d] * csy;
            const double len = std::sqrt(sx * sx + sy * sy);
            stepUx[d] = sx / len;
            stepUy[d] = sy / len;
        }
        const Vec2 wind = windTowardGrid(weather.windFromDeg);
        const double wind20 = weather.wind20ftKmh * KMH_TO_FT_PER_MIN;
        const SurfaceFuelBed fm10 = computeFuelBed(fuelModel10(), weather.moisture);
        constexpr double DEG_TO_RAD = 0.017453292519943295;
        const double minCell = std::min(csx, csy);
        for (int idx = 0; idx < N; ++idx) {
            const Cell& c = grid[idx];
            std::uint8_t fireType = static_cast<std::uint8_t>(FireType::Surface);
            bool burnable = false;
            float head = 0.0f, dur = 0.0f, dir[NUM_DIRS] = {};
            const int code = c.fuelModel;
            const FuelModel* m = code > 0 && code < 256 ? findFuelModel(code) : nullptr;
            const SurfaceFuelBed bed = m ? computeFuelBed(*m, weather.moisture) : SurfaceFuelBed();
            const bool road = ROADS_ARE_FUEL_BREAKS && c.roadClass >= ROAD_CLASS_MIN && c.roadClass <= ROAD_CLASS_MAX;
            if (m && bed.burnable() && !road) {
                const double tanSlope = std::tan(c.slopeDeg * DEG_TO_RAD);
                double waf = weather.windAdjustmentFactor;
                if (waf <= 0.0) {
                    const double crownRatio = c.canopyHeightM > 0.0f
                        ? (c.canopyHeightM - c.canopyBaseHeightM) / c.canopyHeightM : 0.0;
                    waf = windAdjustmentFactor(m->depthFt, c.canopyCoverPct, c.canopyHeightM * FT_PER_M, crownRatio);
                    if (waf != unshelteredWaf(m->depthFt)) ++sheltered;
                }
                SpreadEllipse e = origComputeSpread(bed, wind20 * waf, wind.x, wind.y, tanSlope, c.upslopeX,
                                                    c.upslopeY, weather.windLimit);
                if (weather.crownFire && c.canopyBulkDensity > 0.0f && c.canopyHeightM > 0.0f && c.canopyCoverPct > 0.0f) {
                    const double surfaceIntensity = firelineIntensityKwM(bed, e.headRos);
                    if (surfaceIntensity >= criticalSurfaceIntensity(c.canopyBaseHeightM, weather.foliarMoisturePct)) {
                        SpreadEllipse crown = origComputeSpread(fm10, wind20 * 0.4, wind.x, wind.y, tanSlope,
                                                                c.upslopeX, c.upslopeY, false);
                        crown.headRos *= 3.34;
                        const double crownMPerMin = crown.headRos * FT_PER_MIN_TO_M_PER_MIN;
                        if (crownMPerMin >= criticalActiveRos(c.canopyBulkDensity)) {
                            fireType = static_cast<std::uint8_t>(FireType::Active);
                            if (crown.headRos > e.headRos) e = crown;
                        }
                        else {
                            fireType = static_cast<std::uint8_t>(FireType::Passive);
                        }
                    }
                }
                burnable = true;
                head = static_cast<float>(e.headRos * FT_PER_MIN_TO_M_PER_MIN);
                for (int d = 0; d < NUM_DIRS; ++d)
                    dir[d] = static_cast<float>(e.rosToward(stepUx[d], stepUy[d]) * FT_PER_MIN_TO_M_PER_MIN);
                dur = static_cast<float>(minCell / head + bed.residenceTime);
            }

            // --- Compare bit for bit with what FireSimulation set up ---
            ++cellsCompared;
            if (fireType == static_cast<std::uint8_t>(FireType::Active)) ++active;
            if (fireType == static_cast<std::uint8_t>(FireType::Passive)) ++passive;
            bool same = sim.canBurn(idx) == burnable && sim.grid()[idx].fireType == fireType
                        && sim.headRos(idx) == head;
            if (burnable) {
                same = same && sim.burnDuration(idx) == dur;
                for (int d = 0; d < NUM_DIRS; ++d) same = same && sim.directionalRos(idx, d) == dir[d];
            }
            if (!same) ++differ;
        }
    }
    std::printf("  %ld cells over %zu weathers (%ld active crown, %ld passive, %ld sheltered), %ld differ\n",
                cellsCompared, weathers.size(), active, passive, sheltered, differ);
    CHECK(differ == 0);
    CHECK(active > 0 && passive > 0 && sheltered > 0);   // The crown fire and canopy paths were exercised
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
// Terrain-adjusted wind (wind_solver.h)
// ============================================================================
// These run on 30 m cells with 600 m wind columns (60 x 50 x 20 nodes on the default grid) so each solve takes
// a fraction of a second. The same hill gives the same crest speed-up within ~2 % at the app's 10 m / 100 m.
namespace {
    constexpr float WIND_TEST_CELL_M = 30.0f;

    WindSolverSettings windTestSettings(float alpha = 1.0f) {
        WindSolverSettings s;
        s.columnSizeM = 600.0f;
        s.alpha = alpha;
        return s;
    }

    // A round Gaussian hill, 400 m high with a 2 km standard deviation, in the middle of the grid
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

    // Compass heading (degrees clockwise from north) that a grid-frame vector points to
    double headingDeg(Vec2 v) { return std::atan2(v.x, -v.y) * 180.0 / 3.14159265358979; }
}

// Flat ground with even fuel: the open wind is already mass-conserving, so nothing may change
TEST(terrain_wind_keeps_the_open_wind_on_flat_ground) {
    std::vector<Cell> grid(N);   // Flat, GR2 everywhere
    const WindField f = solveTerrainWind(grid, WIND_TEST_CELL_M, WIND_TEST_CELL_M, 20.0f, 225.0f, windTestSettings());
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
    const WindField f = solveTerrainWind(gaussianHill(), WIND_TEST_CELL_M, WIND_TEST_CELL_M, 20.0f, 270.0f, windTestSettings());
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
    const WindField f = solveTerrainWind(grid, WIND_TEST_CELL_M, WIND_TEST_CELL_M, 20.0f, 225.0f, windTestSettings());
    const double floorHeading = headingDeg(f.toward[at(CX, CY)]);
    std::printf("  (valley floor: %.2f km/h heading %.1f deg)\n", f.speedKmh[at(CX, CY)], floorHeading);
    CHECK(f.converged);
    CHECK(floorHeading > 0.0 && floorHeading < 35.0);   // Turned at least 10 deg towards the valley axis
    CHECK(f.speedKmh[at(CX, CY)] < 0.9f * 20.0f);        // Sheltered floor
}

// Stability: stable air (alpha < 1) is pushed round the hill more than neutral air, unstable air (alpha > 1)
// goes over it and is turned less (WindNinja's convention, alpha = alpha_H / alpha_V)
TEST(terrain_wind_stability_changes_how_much_air_goes_round) {
    const std::vector<Cell> hill = gaussianHill();
    const int d = static_cast<int>(2.0 * 2000.0 / WIND_TEST_CELL_M);
    double turn[3];
    const float alphas[3] = { 0.3f, 1.0f, 3.0f };
    for (int i = 0; i < 3; ++i) {
        const WindField f = solveTerrainWind(hill, WIND_TEST_CELL_M, WIND_TEST_CELL_M, 20.0f, 270.0f, windTestSettings(alphas[i]));
        CHECK(f.converged);
        turn[i] = -f.toward[at(CX - d, CY - d)].y;   // How far north the NW flank wind is turned
    }
    std::printf("  (NW flank turn: stable %.3f, neutral %.3f, unstable %.3f)\n", turn[0], turn[1], turn[2]);
    CHECK(turn[0] > turn[1] && turn[1] > turn[2]);
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
    const FireSimulation plain(grid, w, WIND_TEST_CELL_M, WIND_TEST_CELL_M);
    const FireSimulation withSame(grid, w, WIND_TEST_CELL_M, WIND_TEST_CELL_M, SIM_MINUTES_PER_TICK, &same);
    bool identical = true;
    for (int i = 0; i < N && identical; i += 997)
        for (int dir = 0; dir < NUM_DIRS; ++dir) identical = identical && plain.directionalRos(i, dir) == withSame.directionalRos(i, dir);
    CHECK(identical);

    // A field blowing north: the fire runs north, not east
    WindField north;
    north.toward.assign(N, windTowardGrid(180.0f));   // From the south: blows north (-y)
    north.speedKmh.assign(N, w.wind20ftKmh);
    FireSimulation sim(grid, w, WIND_TEST_CELL_M, WIND_TEST_CELL_M, SIM_MINUTES_PER_TICK, &north);
    sim.ignite(at(250, 250));
    for (int t = 0; t < 1500; ++t) sim.step();
    const int up = reach(sim, 250, 250, 0, -1), east = reach(sim, 250, 250, 1, 0);
    std::printf("  (north %d cells, east %d cells)\n", up, east);
    CHECK(up > 3 * std::max(1, east));
}

// ----------------------------------------------------------------------------
// Wind arrow overlay (wind_overlay.h)
// ----------------------------------------------------------------------------
namespace {
    // A wind field that is the same everywhere: speedKmh, blowing FROM fromDeg
    WindField uniformWind(float speedKmh, float fromDeg) {
        WindField f;
        f.toward.assign(N, windTowardGrid(fromDeg));
        f.speedKmh.assign(N, speedKmh);
        return f;
    }

    // World units per screen pixel for a window as wide as the map (CELL_SIZE pixels per cell at zoom 1)
    float perPixelAt(const Viewport& v) { return v.width() / static_cast<float>(GRID_WIDTH * CELL_SIZE); }
}

// Arrows stay about WIND_ARROW_SPACING_PX apart on screen at every zoom, cover the view, and stay put in map
// space while panning (the lattice is fixed to the map)
TEST(wind_arrows_keep_their_screen_spacing_and_stay_put_when_panning) {
    const WindField f = uniformWind(20.0f, 270.0f);
    for (const float zoom : { 1.0f, 3.0f, 8.0f, 40.0f }) {
        Viewport v(static_cast<float>(GRID_WIDTH * CELL_SIZE), static_cast<float>(GRID_HEIGHT * CELL_SIZE));
        v.zoomAt(v.centreX(), v.centreY(), zoom);
        const float perPixel = perPixelAt(v);
        const float stepPx = windArrowSpacing(perPixel) / perPixel;   // Lattice spacing on screen
        const std::vector<WindArrow> arrows = buildWindArrows(f, v, perPixel, 20.0f, 20.0f);
        // Every arrow's centre lies in the view or within one spacing of it
        const float step = windArrowSpacing(perPixel);
        bool inView = true;
        for (const WindArrow& a : arrows) {
            const float mx = 0.5f * (a.x0 + a.x1), my = 0.5f * (a.y0 + a.y1);
            inView = inView && mx >= v.left() - step && mx <= v.left() + v.width() + step
                            && my >= v.top() - step && my <= v.top() + v.height() + step;
        }
        // About (view / spacing)^2 arrows: the view is covered
        const float across = v.width() / step;
        std::printf("  (zoom %.0f: spacing %.0f px, %zu arrows)\n", zoom, stepPx, arrows.size());
        CHECK(stepPx >= WIND_ARROW_SPACING_PX && stepPx < 2.0f * WIND_ARROW_SPACING_PX + 1e-3f);
        CHECK(inView);
        CHECK(arrows.size() >= static_cast<std::size_t>(across * across * 0.8f));
    }

    // Pan by a third of a spacing: arrows present in both views are at exactly the same map positions
    Viewport v(static_cast<float>(GRID_WIDTH * CELL_SIZE), static_cast<float>(GRID_HEIGHT * CELL_SIZE));
    v.zoomAt(v.centreX(), v.centreY(), 6.0f);
    const float perPixel = perPixelAt(v);
    const std::vector<WindArrow> before = buildWindArrows(f, v, perPixel, 20.0f, 20.0f);
    v.pan(windArrowSpacing(perPixel) / 3.0f, 0.0f);
    const std::vector<WindArrow> after = buildWindArrows(f, v, perPixel, 20.0f, 20.0f);
    int shared = 0, moved = 0;
    for (const WindArrow& a : after) {
        bool found = false;
        for (const WindArrow& b : before) found = found || (a.x0 == b.x0 && a.y0 == b.y0 && a.x1 == b.x1 && a.y1 == b.y1);
        (found ? shared : moved) += 1;
    }
    std::printf("  (after a pan: %d arrows unchanged, %d new at the edge)\n", shared, moved);
    CHECK(shared > 0);
    CHECK(moved <= static_cast<int>(after.size()) / 5);   // Only a new column at the edge
}

// Each arrow points where the wind blows, with its length in proportion to the speed and its colour from the
// speed against the open wind: white = same, blue = slower, magenta = faster
TEST(wind_arrows_show_direction_speed_and_colour) {
    Viewport v(static_cast<float>(GRID_WIDTH * CELL_SIZE), static_cast<float>(GRID_HEIGHT * CELL_SIZE));
    const float perPixel = perPixelAt(v);
    const float step = windArrowSpacing(perPixel);

    // Open wind from the west everywhere: arrows point east (+x), white, 85 % of the spacing long
    const std::vector<WindArrow> east = buildWindArrows(uniformWind(20.0f, 270.0f), v, perPixel, 20.0f, 20.0f);
    CHECK(!east.empty());
    const WindArrow& a = east.front();
    CHECK(a.x1 > a.x0 && std::fabs(a.y1 - a.y0) < 1e-3f);
    CHECK(std::fabs((a.x1 - a.x0) - 0.85f * step) < 1e-3f);
    CHECK(a.hx0 < a.x1 && a.hx1 < a.x1);                     // Head strokes trail back from the tip...
    CHECK((a.hy0 - a.y1) * (a.hy1 - a.y1) < 0.0f);           // ...one either side of the shaft
    CHECK(a.r == 255 && a.g == 255 && a.b == 255);

    // Wind from the south (blows north, -y), half the fastest speed: half the length
    const std::vector<WindArrow> north = buildWindArrows(uniformWind(10.0f, 180.0f), v, perPixel, 20.0f, 20.0f);
    CHECK(north.front().y1 < north.front().y0);
    CHECK(std::fabs((north.front().y0 - north.front().y1) - 0.425f * step) < 1e-3f);

    // Colours: 30 % or more slower = fully blue, 30 % or more faster = fully magenta, in between = partly coloured
    std::uint8_t r, g, b;
    windArrowColour(14.0f, 20.0f, r, g, b);
    CHECK(r == 60 && g == 140 && b == 255);
    windArrowColour(26.0f, 20.0f, r, g, b);
    CHECK(r == 255 && g == 60 && b == 200);
    windArrowColour(23.0f, 20.0f, r, g, b);
    CHECK(r == 255 && g > 60 && g < 255 && b > 200 && b < 255);

    // No field (USE_TERRAIN_WIND off, or nothing solved): no arrows
    CHECK(buildWindArrows(WindField(), v, perPixel, 20.0f, 20.0f).empty());
}

// On the terrain wind over a hill, the arrow at the crest is the longest and magenta (faster than the open
// wind), and the one at the upwind foot is shorter and blue (slower): the overlay shows what the solver found
TEST(wind_arrows_show_the_terrain_wind_over_a_hill) {
    WindSolverSettings s;
    s.columnSizeM = 600.0f;
    const WindField f = solveTerrainWind(gaussianHill(), WIND_TEST_CELL_M, WIND_TEST_CELL_M, 20.0f, 270.0f, s);
    const float maxKmh = *std::max_element(f.speedKmh.begin(), f.speedKmh.end());
    Viewport v(static_cast<float>(GRID_WIDTH * CELL_SIZE), static_cast<float>(GRID_HEIGHT * CELL_SIZE));
    const float perPixel = perPixelAt(v);
    const std::vector<WindArrow> arrows = buildWindArrows(f, v, perPixel, 20.0f, maxKmh);
    // The arrows nearest the crest and nearest the upwind foot (two standard deviations west of it)
    const float cs = static_cast<float>(CELL_SIZE);
    const auto nearest = [&](float x, float y) {
        const WindArrow* best = nullptr;
        float bestD = 1e30f;
        for (const WindArrow& a : arrows) {
            const float dx = 0.5f * (a.x0 + a.x1) - x, dy = 0.5f * (a.y0 + a.y1) - y;
            if (dx * dx + dy * dy < bestD) { bestD = dx * dx + dy * dy; best = &a; }
        }
        return *best;
    };
    const WindArrow crest = nearest(CX * cs, CY * cs);
    const WindArrow foot = nearest((CX - 2.0f * 2000.0f / WIND_TEST_CELL_M) * cs, CY * cs);
    const auto length = [](const WindArrow& a) { return std::hypot(a.x1 - a.x0, a.y1 - a.y0); };
    std::printf("  (crest %.1f km/h colour %d,%d,%d; upwind foot %.1f km/h colour %d,%d,%d)\n",
                crest.speedKmh, crest.r, crest.g, crest.b, foot.speedKmh, foot.r, foot.g, foot.b);
    CHECK(crest.speedKmh > 20.0f && foot.speedKmh < 20.0f);
    CHECK(length(crest) > length(foot));
    CHECK(crest.g < 255 && crest.r == 255);   // Magenta side
    CHECK(foot.r < 255 && foot.b == 255);     // Blue side
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
// between cells are common), with gaps (a ring with a gap plus scattered islands), and on a wide ring like the
// outer zone (with few sites most blocks have a single candidate, which the assignment fills without distance
// tests; with one or two sites, every block or nearly every block does).
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

    // Shape 3: wide ring, like the outer zone round a fire (the only cells the planner splits into regions)
    std::vector<int> outerRing;
    for (int y = 180; y < 820; ++y)
        for (int x = 280; x < 920; ++x) {
            const int r2 = (x - 600) * (x - 600) + (y - 500) * (y - 500);
            if (r2 > 170 * 170 && r2 < 310 * 310) outerRing.push_back(at(x, y));
        }

    for (const auto* cells : { &square, &ring, &outerRing }) {
        for (int sites : { 1, 2, 4, 9, 64 }) {
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


// ----------------------------------------------------------------------------
// Original (pre-2026-11-20) gap fill and perimeter-line trace, kept word for word as references.
// The production versions only skip work that cannot change the answer (gap fill: later passes look only at
// the previous pass's patch cells; trace: each sample's distance is worked out once instead of four times), so
// their output must match these exactly, point for point.
// ----------------------------------------------------------------------------
namespace ref {
    constexpr float RAY_STEP = 0.5f;   // Cells between samples along a lane (two per cell, so no cell is skipped)

    // Distance (cells) from the committed fire at cell (x, y); far away off the grid
    float distAt(const DistanceTransform& edt, int x, int y) {
        if (x < 0 || x >= GRID_WIDTH || y < 0 || y >= GRID_HEIGHT) return 1e9f;
        return static_cast<float>(std::sqrt(edt.squaredDistance(y * GRID_WIDTH + x)));
    }

    // ------------------------------------------------------------------------
    // Downhill direction of the distance field at p (towards the nearest fire), unit length.
    // Central differences two cells either side smooth out the cell-to-cell steps of the field. On the ridge
    // between two fire edges (in a bay) the two sides cancel and the gradient is ill-defined; then `fallback`
    // (the inward normal of the perimeter line, or another sensible direction) is used instead.
    // ------------------------------------------------------------------------
    Vec2 downhill(const DistanceTransform& edt, Vec2 p, Vec2 fallback) {
        const int x = static_cast<int>(std::floor(p.x)), y = static_cast<int>(std::floor(p.y));
        const Vec2 g((distAt(edt, x + 2, y) - distAt(edt, x - 2, y)) * 0.25f,
                     (distAt(edt, x, y + 2) - distAt(edt, x, y - 2)) * 0.25f);
        const float gl = std::sqrt(g.x * g.x + g.y * g.y);
        // An exact distance field has gradient length 1 away from ridges; well below that means a ridge
        // (or the edge of the transform's window, where the far side reads "very far")
        if (gl < 0.5f || gl > 2.0f) {
            const float fl = std::sqrt(fallback.x * fallback.x + fallback.y * fallback.y);
            return fl > 1e-6f ? fallback / fl : Vec2(0.0f, 1.0f);
        }
        return g * (-1.0f / gl);
    }

    // ------------------------------------------------------------------------
    // Inner-band runs along the line p + dir * t for t in [t0, t1], as run end points in increasing t.
    // Sampled every RAY_STEP. When stopAtFire is set the walk also stops once it is clearly past the band:
    // closer to the committed fire than the standoff (the band can't continue on the far side of the fire for
    // this lane), or back out beyond the perimeter line after having been inside.
    // ------------------------------------------------------------------------
    std::vector<Vec2> bandRuns(const DistanceTransform& edt, const std::vector<std::uint8_t>& zone,
                               Vec2 p, Vec2 dir, float t0, float t1, bool stopAtFire) {
        std::vector<Vec2> pts;
        bool inRun = false;
        Vec2 last;
        for (float t = t0; t <= t1 + 1e-4f; t += RAY_STEP) {
            const Vec2 q = p + dir * t;
            const int x = static_cast<int>(std::floor(q.x)), y = static_cast<int>(std::floor(q.y));
            if (x < 0 || x >= GRID_WIDTH || y < 0 || y >= GRID_HEIGHT) {
                if (t > 0.0f) break;   // Ran off the map
                continue;
            }
            if (stopAtFire && t > 0.0f) {
                const float d = distAt(edt, x, y);
                if (d < MIN_FIRE_DISTANCE - 1.0f) break;                    // Past the standoff, at the fire
                if (t > 2.0f && d > PERIMETER_DISTANCE + 1.5f) break;       // Left the band outwards again
            }
            const bool in = zone[y * GRID_WIDTH + x] == 1;
            if (in && !inRun) { pts.push_back(q); inRun = true; }           // A run starts
            if (!in && inRun) { pts.push_back(last); inRun = false; }       // A run ends at the previous sample
            if (in) last = q;
        }
        if (inRun) pts.push_back(last);
        return pts;
    }

    // Distance from point p to segment a-b
    float segDist(Vec2 p, Vec2 a, Vec2 b) {
        const Vec2 ab = b - a, ap = p - a;
        const float l2 = ab.x * ab.x + ab.y * ab.y;
        const float t = l2 > 1e-9f ? std::clamp((ap.x * ab.x + ap.y * ab.y) / l2, 0.0f, 1.0f) : 0.0f;
        return std::sqrt(sqDist(p, a + ab * t));
    }

    // Distance from p to the nearest point of a lane (its runs only, not the gaps between them)
    float laneDist(const SweepLane& lane, Vec2 p) {
        float best = std::numeric_limits<float>::max();
        for (size_t i = 0; i + 1 < lane.pts.size(); i += 2) best = std::min(best, segDist(p, lane.pts[i], lane.pts[i + 1]));
        return best;
    }

    // Marks every cell (inside area, as a local bitmap) whose centre is within LANE_COVER_RADIUS of the lane's
    // runs. Each run is sampled every RAY_STEP and a disc stamped at each sample; a cell within the radius of a
    // sample is within the radius of the run, so a stamped cell is always really covered.
    void stampLane(const SweepLane& lane, const GridRect& area, std::vector<std::uint8_t>& covered) {
        const float r = LANE_COVER_RADIUS, r2 = r * r;
        const int w = area.x1 - area.x0;
        for (size_t i = 0; i + 1 < lane.pts.size(); i += 2) {
            const Vec2 a = lane.pts[i], b = lane.pts[i + 1];
            const float len = std::sqrt(sqDist(a, b));
            const int steps = std::max(1, static_cast<int>(std::ceil(len / RAY_STEP)));
            for (int s = 0; s <= steps; ++s) {
                const Vec2 c = a + (b - a) * (static_cast<float>(s) / steps);
                const int x0 = std::max(area.x0, static_cast<int>(std::floor(c.x - r)));
                const int x1 = std::min(area.x1 - 1, static_cast<int>(std::floor(c.x + r)));
                const int y0 = std::max(area.y0, static_cast<int>(std::floor(c.y - r)));
                const int y1 = std::min(area.y1 - 1, static_cast<int>(std::floor(c.y + r)));
                for (int y = y0; y <= y1; ++y)
                    for (int x = x0; x <= x1; ++x)
                        if (sqDist(Vec2(x + 0.5f, y + 0.5f), c) <= r2) covered[(y - area.y0) * w + (x - area.x0)] = 1;
            }
        }
    }

int fillCoverageGapsOriginal(std::vector<SweepLane>& lanes, const DistanceTransform& perimEdt,
                     const std::vector<std::uint8_t>& zone, const GridRect& area) {
    if (area.empty() || lanes.empty()) return 0;
    const int w = area.x1 - area.x0, h = area.y1 - area.y0;
    std::vector<std::uint8_t> covered(static_cast<size_t>(w) * h, 0);
    for (const SweepLane& l : lanes) stampLane(l, area, covered);

    constexpr int MAX_PASSES = 6;
    int added = 0;
    for (int pass = 0; pass < MAX_PASSES; ++pass) {
        // --- Find patches of unseen inner cells (8-connected flood fill) ---
        std::vector<int> patchOf(static_cast<size_t>(w) * h, -1);
        std::vector<std::vector<int>> patches;   // Local indices of the cells in each patch
        for (int i = 0; i < w * h; ++i) {
            const int gx = area.x0 + i % w, gy = area.y0 + i / w;
            if (covered[i] || patchOf[i] >= 0 || zone[gy * GRID_WIDTH + gx] != 1) continue;
            const int id = static_cast<int>(patches.size());
            patches.emplace_back();
            std::vector<int> stack{ i };
            patchOf[i] = id;
            while (!stack.empty()) {
                const int c = stack.back();
                stack.pop_back();
                patches[id].push_back(c);
                const int cx = c % w, cy = c / w;
                for (int k = 0; k < NUM_DIRS; ++k) {
                    const int nx = cx + DX[k], ny = cy + DY[k];
                    if (nx < 0 || nx >= w || ny < 0 || ny >= h) continue;
                    const int ni = ny * w + nx;
                    if (covered[ni] || patchOf[ni] >= 0) continue;
                    if (zone[(area.y0 + ny) * GRID_WIDTH + area.x0 + nx] != 1) continue;
                    patchOf[ni] = id;
                    stack.push_back(ni);
                }
            }
        }
        if (patches.empty()) break;   // Everything is seen

        // --- One short lane through each patch ---
        const size_t mainCount = lanes.size();
        std::vector<SweepLane> extra;
        for (const std::vector<int>& patch : patches) {
            // Patch centre: the patch cell nearest the mean of its cells (the mean itself can lie outside a
            // curved patch)
            Vec2 mean(0.0f, 0.0f);
            for (int c : patch) mean += Vec2(area.x0 + c % w + 0.5f, area.y0 + c / w + 0.5f);
            mean = mean / static_cast<float>(patch.size());
            Vec2 centre = mean;
            float bd = std::numeric_limits<float>::max();
            for (int c : patch) {
                const Vec2 p(area.x0 + c % w + 0.5f, area.y0 + c / w + 0.5f);
                if (sqDist(p, mean) < bd) { bd = sqDist(p, mean); centre = p; }
            }

            // Same direction as the main lanes: down the distance-field gradient at the centre. Fallback on a
            // ridge: towards the nearest main lane's inner end direction (any lane direction will do there).
            const SweepLane& nearLane = *std::min_element(lanes.begin(), lanes.begin() + mainCount,
                [&](const SweepLane& a, const SweepLane& b) { return laneDist(a, centre) < laneDist(b, centre); });
            const Vec2 laneDir = nearLane.pts.back() - nearLane.pts.front();
            const Vec2 dir = downhill(perimEdt, centre, laneDir.x * laneDir.x + laneDir.y * laneDir.y > 1e-6f
                                                         ? laneDir : Vec2(0.0f, 1.0f));

            // How far along that line the patch reaches (cells close enough to the line to be seen from it)
            float tMin = 0.0f, tMax = 0.0f;
            for (int c : patch) {
                const Vec2 rel = Vec2(area.x0 + c % w + 0.5f, area.y0 + c / w + 0.5f) - centre;
                const float t = rel.x * dir.x + rel.y * dir.y;
                const float perp = std::fabs(rel.x * dir.y - rel.y * dir.x);
                if (perp > LANE_COVER_RADIUS) continue;
                tMin = std::min(tMin, t);
                tMax = std::max(tMax, t);
            }

            SweepLane lane;
            lane.pts = bandRuns(perimEdt, zone, centre, dir, tMin - RAY_STEP, tMax + RAY_STEP, false);
            if (lane.pts.empty()) lane.pts = { centre, centre };   // Centre is an inner cell, so can't happen
            lane.seed = lane.pts.front();
            lane.loop = nearLane.loop;
            lane.gapFill = true;
            extra.push_back(std::move(lane));
        }

        // --- Insert each new lane after the main lane it sits closest to (and after earlier gap lanes there) ---
        for (SweepLane& g : extra) {
            stampLane(g, area, covered);
            const Vec2 mid = (g.pts.front() + g.pts.back()) * 0.5f;
            size_t parent = 0;
            float pd = std::numeric_limits<float>::max();
            for (size_t i = 0; i < lanes.size(); ++i) {
                if (lanes[i].gapFill || lanes[i].loop != g.loop) continue;
                const float d = laneDist(lanes[i], mid);
                if (d < pd) { pd = d; parent = i; }
            }
            size_t at = parent + 1;
            while (at < lanes.size() && lanes[at].gapFill) ++at;
            lanes.insert(lanes.begin() + static_cast<std::ptrdiff_t>(at), std::move(g));
            ++added;
        }
    }
    return added;
}

    // Corner offsets of a square, in walking order: top-left, top-right, bottom-right, bottom-left
    constexpr int SQ_DX[4] = { 0, 1, 1, 0 };
    constexpr int SQ_DY[4] = { 0, 0, 1, 1 };

    // Id of the square side from corner c to corner c+1 (mod 4) of the square whose top-left sample is (x, y).
    // Horizontal sides and vertical sides get even and odd ids, keyed by the side's top / left sample, so the
    // side shared by two squares gets the same id from both. Samples run from -1 to GRID_WIDTH (one ring of
    // padding round the grid), hence the +1 offsets.
    std::int64_t sideId(int x, int y, int c) {
        const std::int64_t W = GRID_WIDTH + 2;
        auto key = [&](int sx, int sy, int vertical) {
            return ((static_cast<std::int64_t>(sy) + 1) * W + (sx + 1)) * 2 + vertical;
        };
        switch (c) {
        case 0:  return key(x, y, 0);         // Top side
        case 1:  return key(x + 1, y, 1);     // Right side
        case 2:  return key(x, y + 1, 0);     // Bottom side
        default: return key(x, y, 1);         // Left side
        }
    }

    // One traced piece of line inside a square, from one side crossing to another
    struct Piece {
        Vec2 from, to;
        std::int64_t toSide;   // Side the piece ends on, i.e. where the next piece starts
        bool used = false;
    };

std::vector<std::vector<Vec2>> tracePerimeterLineOriginal(const DistanceTransform& perimEdt, float level) {
    // Only the transform's window holds finite distances; the level line lies well inside it, so squares
    // from one sample outside the window to its far edge are enough
    const GridRect win = perimEdt.window();
    std::vector<std::vector<Vec2>> loops;
    if (win.empty()) return loops;

    // Level-line function at sample (x, y): negative on the fire side of the line (nearer than `level`).
    // Samples off the grid count as far from the fire, so a loop that would leave the map closes along its edge.
    auto f = [&](int x, int y) -> double {
        if (x < 0 || x >= GRID_WIDTH || y < 0 || y >= GRID_HEIGHT) return 1e9;
        return std::sqrt(perimEdt.squaredDistance(y * GRID_WIDTH + x)) - level;
    };

    // Pieces keyed by the side they start on (each side crossing starts exactly one piece)
    std::unordered_map<std::int64_t, Piece> pieces;
    std::vector<std::int64_t> order;   // Start sides in scan order, so tracing is deterministic

    for (int y = win.y0 - 1; y < win.y1; ++y) {
        for (int x = win.x0 - 1; x < win.x1; ++x) {
            double v[4];
            bool in[4];
            int nIn = 0;
            for (int c = 0; c < 4; ++c) {
                v[c] = f(x + SQ_DX[c], y + SQ_DY[c]);
                in[c] = v[c] < 0.0;
                nIn += in[c];
            }
            if (nIn == 0 || nIn == 4) continue;   // The line doesn't pass through this square

            // Crossing point on side c (corner c -> c+1), and whether the walk leaves the fire side there
            Vec2 at[4];
            bool leaving[4], crosses[4];
            for (int c = 0; c < 4; ++c) {
                const int d = (c + 1) % 4;
                crosses[c] = in[c] != in[d];
                leaving[c] = in[c] && !in[d];
                if (!crosses[c]) continue;
                const double t = v[c] / (v[c] - v[d]);   // Where the linear interpolation reaches 0
                const Vec2 a(x + SQ_DX[c] + 0.5f, y + SQ_DY[c] + 0.5f), b(x + SQ_DX[d] + 0.5f, y + SQ_DY[d] + 0.5f);
                at[c] = a + (b - a) * static_cast<float>(t);
            }

            // Pair each leaving crossing with an entering one. With two crossings there is only one way. With
            // four (a saddle: opposite corners on the same side) the centre value decides: if the centre is on
            // the fire side, the two fire corners are joined, so each leaving crossing pairs with the NEXT
            // entering crossing round the square (cutting off the far corners); otherwise with the previous one.
            const bool saddle = nIn == 2 && in[0] == in[2];
            const bool centreIn = (v[0] + v[1] + v[2] + v[3]) * 0.25 < 0.0;
            for (int c = 0; c < 4; ++c) {
                if (!crosses[c] || !leaving[c]) continue;
                int e = -1;
                if (saddle && !centreIn) {
                    for (int k = 1; k < 4 && e < 0; ++k) { const int j = (c + 4 - k) % 4; if (crosses[j] && !leaving[j]) e = j; }
                }
                else {
                    for (int k = 1; k < 4 && e < 0; ++k) { const int j = (c + k) % 4; if (crosses[j] && !leaving[j]) e = j; }
                }
                if (e < 0) continue;   // Cannot happen with a consistent sign pattern
                const std::int64_t start = sideId(x, y, c);
                pieces[start] = Piece{ at[c], at[e], sideId(x, y, e) };
                order.push_back(start);
            }
        }
    }

    // Chain pieces into loops: each piece ends on the side where the next one starts
    for (std::int64_t first : order) {
        auto it = pieces.find(first);
        if (it == pieces.end() || it->second.used) continue;
        std::vector<Vec2> loop;
        while (it != pieces.end() && !it->second.used) {
            it->second.used = true;
            loop.push_back(it->second.from);
            it = pieces.find(it->second.toSide);
        }
        if (loop.size() >= 3) loops.push_back(std::move(loop));
    }

    // Make every loop run anticlockwise on screen (negative area in y-down grid coordinates). The piece rule
    // gives every loop the same sense round the fire, so if the first enclosing loop is the wrong way round,
    // they all are: test the loop with the largest area (always one that encloses fire) and flip all if needed.
    if (!loops.empty()) {
        size_t biggest = 0;
        for (size_t i = 1; i < loops.size(); ++i) {
            if (std::fabs(loopSignedArea(loops[i])) > std::fabs(loopSignedArea(loops[biggest]))) biggest = i;
        }
        if (loopSignedArea(loops[biggest]) > 0.0) {
            for (auto& l : loops) std::reverse(l.begin(), l.end());
        }
    }

    // Start each loop at its topmost point (then leftmost), a stable place for lane numbering to begin
    for (auto& l : loops) {
        const auto top = std::min_element(l.begin(), l.end(), [](Vec2 a, Vec2 b) {
            return a.y < b.y || (a.y == b.y && a.x < b.x);
        });
        std::rotate(l.begin(), top, l.end());
    }
    return loops;
}
}   // namespace ref

namespace {
    // True if two lane lists are identical: same lanes in the same order, every point bit for bit
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

    // Fires for the exactness tests below: the bay (gap fill needs several passes), a fire with two bays, a fire
    // on the grid's west edge (its perimeter line runs off the map and closes along the edge), a random
    // blobby fire (many small bays and patches) and two separate fires (two loops)
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
}

// Lane coverage is now worked out per cell (only the samples that can reach each cell are tested). It must mark
// exactly the cells the original disc-per-sample stamp (ref::stampLane, kept above) marks. Random lanes at every
// angle and length, with several runs, axis-aligned and zero-length runs, clipped by the area's edges.
TEST(lane_stamp_matches_disc_stamps) {
    std::mt19937 rng(11);
    std::uniform_real_distribution<float> pos(300.0f, 700.0f), len(0.0f, 300.0f), ang(0.0f, 6.2831853f), gap(0.0f, 40.0f);
    const GridRect area{ 380, 330, 760, 640 };   // Smaller than the spread of the lanes, so many are clipped
    const size_t cells = static_cast<size_t>(area.x1 - area.x0) * (area.y1 - area.y0);
    int lanes = 0, mismatched = 0;
    long long marked = 0;
    for (int k = 0; k < 2000; ++k) {
        SweepLane lane;
        Vec2 p(pos(rng), pos(rng));
        // Direction: every 10th lane axis-aligned (a special case of the row interval), the rest at random
        const float th = k % 10 == 0 ? 1.5707963f * (k / 10 % 4) : ang(rng);
        const Vec2 dir(std::cos(th), std::sin(th));
        const int runs = 1 + k % 3;
        for (int r = 0; r < runs; ++r) {
            // Every 7th run has zero length, every 13th a tiny one (just above the zero-length cut-off)
            const float l = k % 7 == 0 ? 0.0f : (k % 13 == 0 ? 2e-3f : len(rng) / runs);
            lane.pts.push_back(p);
            p = p + dir * l;
            lane.pts.push_back(p);
            p = p + dir * gap(rng);
        }
        std::vector<std::uint8_t> fast(cells, 0), original(cells, 0);
        stampLaneCoverage(lane, area, fast);
        ref::stampLane(lane, area, original);
        mismatched += fast != original;
        for (std::uint8_t v : original) marked += v;
        ++lanes;
    }
    std::printf("  (%d random lanes, %lld cells marked, %d lanes differ)\n", lanes, marked, mismatched);
    CHECK(mismatched == 0);
    CHECK(marked > 0);
}

// The stamp now tests a cell only against the samples within h = sqrt((r + 0.1)^2 - across^2) of its projection
// onto the run (2026-11-23). The cells where skipping a sample could go wrong are the ones right on the swath's
// edge, where the nearest sample is almost exactly LANE_COVER_RADIUS away. Lanes on whole and half-cell
// coordinates put cell centres EXACTLY r = 6.5 cells across from the run (and exactly on samples along it), and
// lanes nudged by a few float ulps put them a hair inside or outside; long lanes far out on a 2000-cell grid
// have the largest rounding. Every stamp must match the disc stamps bit for bit.
TEST(lane_stamp_matches_disc_stamps_on_the_edge) {
    std::mt19937 rng(23);
    std::uniform_int_distribution<int> cell(20, 1950);
    std::uniform_int_distribution<int> lenCells(0, 400), nudge(-3, 3), half(0, 1), axis(0, 5);
    const GridRect area{ 0, 0, 2000, 2000 };   // Whole 2000 x 2000 grid, coordinates up to 2000 (most rounding)
    const size_t cells = static_cast<size_t>(area.x1 - area.x0) * (area.y1 - area.y0);
    std::vector<std::uint8_t> fast(cells), original(cells);
    int lanes = 0, mismatched = 0;
    long long marked = 0;
    // Moves a float by k units in the last place (k may be negative), to sit just either side of the edge
    auto ulps = [](float v, int k) {
        for (; k > 0; --k) v = std::nextafter(v, 1e9f);
        for (; k < 0; ++k) v = std::nextafter(v, -1e9f);
        return v;
    };
    for (int k = 0; k < 600; ++k) {
        // Start on a whole or half cell, nudged a few ulps; run along an axis, a diagonal or a 1:2 slope
        const float sx = ulps(cell(rng) + 0.5f * half(rng), nudge(rng));
        const float sy = ulps(cell(rng) + 0.5f * half(rng), nudge(rng));
        static const int DIRS[6][2] = { {1, 0}, {0, 1}, {-1, 0}, {1, 1}, {1, -2}, {2, 1} };
        const int* d = DIRS[axis(rng)];
        const float l = static_cast<float>(lenCells(rng));
        const float norm = std::sqrt(static_cast<float>(d[0] * d[0] + d[1] * d[1]));
        const Vec2 a(sx, sy);
        Vec2 b(std::clamp(sx + d[0] * l / norm, 1.0f, 1999.0f), std::clamp(sy + d[1] * l / norm, 1.0f, 1999.0f));
        b = Vec2(ulps(b.x, nudge(rng)), ulps(b.y, nudge(rng)));
        SweepLane lane;
        lane.pts = { a, b };
        std::fill(fast.begin(), fast.end(), 0);
        std::fill(original.begin(), original.end(), 0);
        stampLaneCoverage(lane, area, fast);
        ref::stampLane(lane, area, original);
        mismatched += fast != original;
        for (std::uint8_t v : original) marked += v;
        ++lanes;
    }
    std::printf("  (%d edge lanes, %lld cells marked, %d lanes differ)\n", lanes, marked, mismatched);
    CHECK(mismatched == 0);
    CHECK(marked > 0);
}

// Gap fill now only re-floods the previous pass's patches (lanes only remove gaps, never add them). The lanes
// it adds must be exactly the ones the original whole-window scan added, in the same places and order.
TEST(gap_fill_matches_full_window_scan) {
    int totalAdded = 0;
    for (const auto& fire : exactnessFires()) {
        SweepSetup s(fire);
        // Rebuild the gap fill's input the way planMission does, then run both versions on copies of it
        const std::vector<SweepLane> lanes = layLanes(tracePerimeterLine(s.perim.distance()), s.perim.distance(),
                                                      s.plan.zone);
        std::vector<SweepLane> fast = lanes, original = lanes;
        const int addedFast = fillCoverageGaps(fast, s.perim.distance(), s.plan.zone, s.plan.bounds);
        const int addedOriginal = ref::fillCoverageGapsOriginal(original, s.perim.distance(), s.plan.zone, s.plan.bounds);
        CHECK(addedFast == addedOriginal);
        CHECK(sameLanes(fast, original));
        totalAdded += addedOriginal;
    }
    std::printf("  (%d gap-fill lanes added over the test fires, identical in both versions)\n", totalAdded);
    CHECK(totalAdded > 0);   // The fires really exercise the gap fill
}

// Original (pre-2026-11-22) ownership flood, kept word for word: seed every cell a lane run passes over, then
// breadth-first through the band, bounds-testing every neighbour. The last stage of buildOwnerMap (band cells cut
// off from every lane go to the nearest lane) is unchanged and left out here, so those cells stay -1.
namespace ref {
    void ownerSeedAndFlood(const std::vector<SweepLane>& lanes, const std::vector<std::uint8_t>& zone,
                           const GridRect& area, std::vector<int>& owner) {
        constexpr float RAY_STEP = 0.5f;
        std::vector<int> queue;
        for (const SweepLane& l : lanes) {
            for (size_t i = 0; i + 1 < l.pts.size(); i += 2) {
                const Vec2 a = l.pts[i], b = l.pts[i + 1];
                const int steps = std::max(1, static_cast<int>(std::ceil(std::sqrt(sqDist(a, b)) / RAY_STEP)));
                for (int s = 0; s <= steps; ++s) {
                    const Vec2 p = a + (b - a) * (static_cast<float>(s) / steps);
                    const int x = static_cast<int>(std::floor(p.x)), y = static_cast<int>(std::floor(p.y));
                    if (x < area.x0 || x >= area.x1 || y < area.y0 || y >= area.y1) continue;
                    const int idx = y * GRID_WIDTH + x;
                    if (zone[idx] != 1 || owner[idx] >= 0) continue;
                    owner[idx] = l.sector;
                    queue.push_back(idx);
                }
            }
        }
        for (size_t head = 0; head < queue.size(); ++head) {
            const int idx = queue[head];
            const int x = idx % GRID_WIDTH, y = idx / GRID_WIDTH;
            for (int k = 0; k < NUM_DIRS; ++k) {
                const int nx = x + DX[k], ny = y + DY[k];
                if (nx < area.x0 || nx >= area.x1 || ny < area.y0 || ny >= area.y1) continue;
                const int ni = ny * GRID_WIDTH + nx;
                if (zone[ni] != 1 || owner[ni] >= 0) continue;
                owner[ni] = owner[idx];
                queue.push_back(ni);
            }
        }
    }
}   // namespace ref

// buildOwnerMap's breadth-first flood now skips the bounds tests for cells inside the area and reserves its queue.
// Every cell the flood reaches must get exactly the sector the original flood gave it (same lane reaches it first),
// and the cells it cannot reach must still be given a sector by the last stage. Run on the test fires with the
// planner's own window and with smaller areas cut out of it, so the area's edge runs through the band (many edge
// cells), plus an island of band cells no lane reaches.
TEST(owner_map_matches_original_flood) {
    long long compared = 0, cutOff = 0, differ = 0, unowned = 0;
    for (const auto& fire : exactnessFires()) {
        SweepSetup s(fire);
        const GridRect w = s.plan.bounds;
        const int mx = (w.x0 + w.x1) / 2, my = (w.y0 + w.y1) / 2;
        // An island of band cells in the window's corner, well away from the band and every lane, so the flood
        // cannot reach it and the last stage has to place it
        std::vector<std::uint8_t> zone = s.plan.zone;
        for (int y = w.y0 + 2; y < w.y0 + 12; ++y)
            for (int x = w.x0 + 2; x < w.x0 + 12; ++x) zone[at(x, y)] = 1;
        for (const GridRect area : { w, GridRect{ w.x0, w.y0, mx, my }, GridRect{ mx - 37, w.y0, w.x1, my + 23 },
                                     GridRect{ w.x0 + 5, my - 60, w.x1 - 5, my + 60 } }) {
            std::vector<int> fast(N, -1), original(N, -1);
            buildOwnerMap(s.plan.lanes, zone, area, fast);
            ref::ownerSeedAndFlood(s.plan.lanes, zone, area, original);
            for (int y = area.y0; y < area.y1; ++y) {
                for (int x = area.x0; x < area.x1; ++x) {
                    const int i = at(x, y);
                    if (zone[i] != 1) { differ += fast[i] != -1; continue; }
                    if (original[i] >= 0) { differ += fast[i] != original[i]; ++compared; }
                    else { ++cutOff; unowned += fast[i] < 0 && !s.plan.lanes.empty(); }
                }
            }
            // Nothing outside the area is written
            for (int i = 0; i < N; ++i) {
                const int x = i % GRID_WIDTH, y = i / GRID_WIDTH;
                if (x < area.x0 || x >= area.x1 || y < area.y0 || y >= area.y1) differ += fast[i] != -1;
            }
        }
    }
    std::printf("  (%lld flooded cells compared, %lld cut off from every lane; %lld differ, %lld left unowned)\n",
                compared, cutOff, differ, unowned);
    CHECK(differ == 0);
    CHECK(unowned == 0);
    CHECK(compared > 100000);
    CHECK(cutOff > 0);   // The cut-down areas really leave band pieces the flood cannot reach
}

// buildOwnerMap can be told how many band cells the area holds (the planner counts them while classifying) and
// then skips its final scan for cut-off cells when the flood has already reached them all. The map must be
// identical with and without the count, both when nothing is cut off (scan skipped) and when an island of band
// cells is cut off from every lane (scan still runs).
TEST(owner_map_band_count_skip_is_exact) {
    int skippedCases = 0, scannedCases = 0;
    long long differ = 0;
    for (const auto& fire : exactnessFires()) {
        SweepSetup s(fire);
        const GridRect w = s.plan.bounds;
        std::vector<std::uint8_t> island = s.plan.zone;   // As owner_map_matches_original_flood: a far island
        for (int y = w.y0 + 2; y < w.y0 + 12; ++y)
            for (int x = w.x0 + 2; x < w.x0 + 12; ++x) island[at(x, y)] = 1;
        for (const std::vector<std::uint8_t>* zone : { &s.plan.zone, &island }) {
            long long band = 0;
            for (int y = w.y0; y < w.y1; ++y)
                for (int x = w.x0; x < w.x1; ++x) band += (*zone)[at(x, y)] == 1;
            const std::vector<int> plain = buildOwnerMap(s.plan.lanes, *zone, w);          // Always scans
            const std::vector<int> counted = buildOwnerMap(s.plan.lanes, *zone, w, band);  // May skip the scan
            for (int i = 0; i < N; ++i) differ += plain[i] != counted[i];
            // Every band cell has an owner either way (so a skipped scan really had nothing left to do)
            long long unowned = 0;
            for (int y = w.y0; y < w.y1; ++y)
                for (int x = w.x0; x < w.x1; ++x) unowned += (*zone)[at(x, y)] == 1 && counted[at(x, y)] < 0;
            CHECK(unowned == 0 || s.plan.lanes.empty());
            (zone == &island ? scannedCases : skippedCases) += 1;
        }
    }
    std::printf("  (%d plain and %d cut-off cases; %lld cells differ)\n", skippedCases, scannedCases, differ);
    CHECK(differ == 0);
    CHECK(skippedCases > 0 && scannedCases > 0);
}

// The planner now builds the ownership map with buildOwnerMapUsingZone, which marks claimed band cells in the
// zone map itself (one byte read per neighbour instead of the zone and owner arrays) and puts the marks back.
// Checked on the test fires, with the planner's window and cut-down areas (many edge cells, pieces cut off from
// every lane) and a far island of band cells:
//  - every cell the flood reaches gets exactly the sector the original flood (ref::ownerSeedAndFlood) gave it;
//  - cut-off band cells still get an owner from the final scan, and nothing outside the area is written;
//  - the zone map is exactly as it was on entry, with and without the band count (which skips the final scan).
TEST(owner_map_using_zone_matches_and_restores_zone) {
    long long compared = 0, differ = 0, unowned = 0, zoneChanged = 0;
    int cases = 0;
    for (const auto& fire : exactnessFires()) {
        SweepSetup s(fire);
        const GridRect w = s.plan.bounds;
        const int mx = (w.x0 + w.x1) / 2, my = (w.y0 + w.y1) / 2;
        std::vector<std::uint8_t> island = s.plan.zone;   // As owner_map_matches_original_flood: a far island
        for (int y = w.y0 + 2; y < w.y0 + 12; ++y)
            for (int x = w.x0 + 2; x < w.x0 + 12; ++x) island[at(x, y)] = 1;
        for (const std::vector<std::uint8_t>* zoneIn : { &s.plan.zone, &island }) {
            for (const GridRect area : { w, GridRect{ w.x0, w.y0, mx, my }, GridRect{ mx - 37, w.y0, w.x1, my + 23 } }) {
                // Band cells in the area, for the form that may skip the final scan
                long long band = 0;
                for (int y = area.y0; y < area.y1; ++y)
                    for (int x = area.x0; x < area.x1; ++x) band += (*zoneIn)[at(x, y)] == 1;
                std::vector<int> original(N, -1);
                ref::ownerSeedAndFlood(s.plan.lanes, *zoneIn, area, original);
                for (const long long count : { -1LL, band }) {
                    std::vector<std::uint8_t> zone = *zoneIn;   // Marked and restored by the call
                    std::vector<int> fast(N, -1);
                    buildOwnerMapUsingZone(s.plan.lanes, zone, area, fast, count);
                    zoneChanged += zone != *zoneIn;
                    ++cases;
                    for (int i = 0; i < N; ++i) {
                        const int x = i % GRID_WIDTH, y = i / GRID_WIDTH;
                        const bool inside = x >= area.x0 && x < area.x1 && y >= area.y0 && y < area.y1;
                        if (!inside || (*zoneIn)[i] != 1) { differ += fast[i] != -1; continue; }
                        if (original[i] >= 0) { differ += fast[i] != original[i]; ++compared; }
                        else unowned += fast[i] < 0 && !s.plan.lanes.empty();
                    }
                }
            }
        }
    }
    std::printf("  (%d cases, %lld flooded cells compared; %lld differ, %lld left unowned, zone changed in %lld)\n",
                cases, compared, differ, unowned, zoneChanged);
    CHECK(differ == 0);
    CHECK(unowned == 0);
    CHECK(zoneChanged == 0);
    CHECK(compared > 100000);
}

// Original (pre-2026-11-25) sector matching, kept word for word: it counted shared cells by scanning every cell of
// the plan window and searched the fleet's ids for each one. The planner now counts over its list of band cells
// with a sector-to-drone table worked out once.
namespace ref {
    std::vector<int> matchSectorsOriginal(const MissionPlan& plan, const std::vector<std::vector<Vec2>>& routes,
                                          const SweepFleet& fleet, const GridRect& area) {
        const int S = static_cast<int>(routes.size()), D = static_cast<int>(fleet.ids.size());
        std::vector<int> result(S, -1);
        std::vector<char> droneUsed(D, 0);
        struct Pair { double score; int sector, drone; };
        std::vector<Pair> pairs;
        const MissionPlan* prev = fleet.previous;
        if (prev && !prev->owner.empty() && !prev->sectorDrone.empty()) {
            std::vector<long long> shared(static_cast<size_t>(S) * D, 0);
            for (int y = area.y0; y < area.y1; ++y) {
                for (int i = y * GRID_WIDTH + area.x0, end = y * GRID_WIDTH + area.x1; i < end; ++i) {
                    const int s = plan.owner[i], ps = prev->owner[i];
                    if (s < 0 || ps < 0 || ps >= static_cast<int>(prev->sectorDrone.size())) continue;
                    const int id = prev->sectorDrone[ps];
                    for (int d = 0; d < D; ++d) if (fleet.ids[d] == id) { ++shared[static_cast<size_t>(s) * D + d]; break; }
                }
            }
            for (int s = 0; s < S; ++s)
                for (int d = 0; d < D; ++d)
                    if (shared[static_cast<size_t>(s) * D + d] > 0) pairs.push_back({ static_cast<double>(shared[static_cast<size_t>(s) * D + d]), s, d });
            std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) {
                return a.score > b.score || (a.score == b.score && (a.sector < b.sector || (a.sector == b.sector && a.drone < b.drone)));
            });
            for (const Pair& p : pairs) {
                if (result[p.sector] >= 0 || droneUsed[p.drone]) continue;
                result[p.sector] = fleet.ids[p.drone];
                droneUsed[p.drone] = 1;
            }
        }
        pairs.clear();
        for (int s = 0; s < S; ++s) {
            if (result[s] >= 0 || routes[s].empty()) continue;
            for (int d = 0; d < D; ++d) {
                if (droneUsed[d]) continue;
                float best = std::numeric_limits<float>::max();
                for (const Vec2& p : routes[s]) best = std::min(best, sqDist(p, fleet.positions[d]));
                pairs.push_back({ best, s, d });
            }
        }
        std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) {
            return a.score < b.score || (a.score == b.score && (a.sector < b.sector || (a.sector == b.sector && a.drone < b.drone)));
        });
        for (const Pair& p : pairs) {
            if (result[p.sector] >= 0 || droneUsed[p.drone]) continue;
            result[p.sector] = fleet.ids[p.drone];
            droneUsed[p.drone] = 1;
        }
        for (int s = 0; s < S; ++s) {
            if (result[s] >= 0) continue;
            for (int d = 0; d < D; ++d) if (!droneUsed[d]) { result[s] = fleet.ids[d]; droneUsed[d] = 1; break; }
        }
        return result;
    }
}   // namespace ref

// Sector matching now counts shared land over the band cells only, with each previous sector's drone looked up
// once. Over a growing, moving fire, each plan's sectorDrone must equal what the original window scan gives for the
// same plan and fleet. Three kinds of fleet are tried at every step: the previous plan's drones in order; the same
// drones reversed (so fleet index and id differ); and one drone retired plus one the previous plan never saw (so
// some previous sectors have no drone, and one drone has no previous land).
TEST(sector_matching_over_band_matches_window_scan) {
    const std::vector<std::vector<std::uint8_t>> fires = {
        discFire(CX, CY, 40), discFire(CX, CY, 55), bayFire(), discFire(CX + 30, CY - 20, 70), discFire(CX + 30, CY - 20, 90) };
    int compared = 0, differ = 0, matchedByLand = 0;
    MissionPlan prev;
    int nextId = 100;
    for (const auto& fire : fires) {
        DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
        live.compute(fire);
        Perimeter perim(fire);
        SafePlanner sp(live);
        MissionPlan kept;   // The plan the next fire's fleets follow on from
        for (int variant = 0; variant < 3; ++variant) {
            SweepFleet fleet;
            if (prev.sectorDrone.empty()) {
                for (int i = 0; i < SWEEP_DRONES; ++i) { fleet.ids.push_back(i); fleet.positions.push_back(Vec2(0.5f, 0.5f)); }
            }
            else {
                for (size_t k = 0; k < prev.sectorDrone.size(); ++k) {
                    fleet.ids.push_back(prev.sectorDrone[k]);
                    fleet.positions.push_back(prev.sweepRoutes[k].empty() ? Vec2(0.5f, 0.5f) : prev.sweepRoutes[k].front());
                }
                if (variant == 1) {
                    std::reverse(fleet.ids.begin(), fleet.ids.end());
                    std::reverse(fleet.positions.begin(), fleet.positions.end());
                }
                if (variant == 2 && fleet.ids.size() > 1) {
                    fleet.ids.erase(fleet.ids.begin());             // Retired: its old sector has no drone now
                    fleet.positions.erase(fleet.positions.begin());
                    fleet.ids.push_back(999);                       // New: no land in the previous plan
                    fleet.positions.push_back(Vec2(CX + 0.5f, CY + 0.5f));
                }
                fleet.previous = &prev;
            }
            std::mt19937 rng(2024);
            MissionPlan plan = planMission(live, perim.distance(), sp, rng, fleet);
            const std::vector<int> original = ref::matchSectorsOriginal(plan, plan.sweepRoutes, fleet, plan.bounds);
            differ += plan.sectorDrone != original;
            ++compared;
            // Sectors given a drone already in the fleet (by shared land or distance), so matching really ran
            if (fleet.previous)
                for (int id : plan.sectorDrone) matchedByLand += std::find(fleet.ids.begin(), fleet.ids.end(), id) != fleet.ids.end();
            if (variant == 0) kept = std::move(plan);
        }
        for (int& id : kept.sectorDrone) if (id < 0) id = nextId++;   // As the controller does: new drones get ids
        prev = std::move(kept);
    }
    std::printf("  (%d plans compared, %d differ; %d sectors given a fleet drone)\n", compared, differ, matchedByLand);
    CHECK(differ == 0);
    CHECK(compared == 15);
    CHECK(matchedByLand > 0);
}

// The perimeter-line trace now works out each sample once and slides two row buffers down the window. The
// loops must be exactly those of the original, which worked out all four corners of every square.
TEST(perimeter_trace_matches_per_corner_evaluation) {
    int loopsChecked = 0;
    for (const auto& fire : exactnessFires()) {
        Perimeter perim(fire);
        // The planner's level, the survey standoff, and an odd fractional level
        for (float level : { static_cast<float>(PERIMETER_DISTANCE), static_cast<float>(MIN_FIRE_DISTANCE), 57.3f }) {
            const auto fast = tracePerimeterLine(perim.distance(), level);
            const auto original = ref::tracePerimeterLineOriginal(perim.distance(), level);
            CHECK(fast.size() == original.size());
            if (fast.size() != original.size()) continue;
            for (size_t l = 0; l < fast.size(); ++l) {
                bool same = fast[l].size() == original[l].size();
                for (size_t k = 0; same && k < fast[l].size(); ++k)
                    same = fast[l][k].x == original[l][k].x && fast[l][k].y == original[l][k].y;
                CHECK(same);
                ++loopsChecked;
            }
        }
    }
    std::printf("  (%d loops identical point for point)\n", loopsChecked);
    CHECK(loopsChecked > 0);
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

// The sector split now reads its laps from a SectorLapTable, which remembers them and finds each greedy sector's
// last lane by binary search over the running maximum of the laps. Both must match the original code bit for
// bit: every lap(first, last) must equal estimateSectorLength, and greedyLast must return exactly the lane the
// original lane-by-lane fill stopped at (grow a SectorLapBuilder while lapWithNext() <= limit). Limits are
// chosen to sit exactly on lap values and one float step either side, where a rounding slip would show, as well
// as at random; questions come in random order so rows are grown lazily in every pattern; and the table is reset
// and reused after the lanes change, as the planner does after joining loops.
TEST(sector_lap_table_matches_lane_by_lane_fill) {
    // Reference: the greedy fill as it was in splitGroup before the table
    auto greedyRef = [](const std::vector<SweepLane>& lanes, int first, int bound, float limit) {
        SectorLapBuilder b(lanes, first);
        while (b.last() + 1 <= bound && b.lapWithNext() <= limit) b.addNext();
        return b.last();
    };

    std::mt19937 rng(11);
    std::uniform_real_distribution<float> coord(0.0f, 1500.0f);
    auto randomLanes = [&](int n) {
        std::vector<SweepLane> lanes(n);
        for (SweepLane& l : lanes) {
            const int pts = 2 * (1 + static_cast<int>(rng() % 3));
            for (int k = 0; k < pts; ++k) l.pts.push_back(Vec2(coord(rng), coord(rng)));
            if (rng() % 5 == 0) l.pts[1] = l.pts[0];   // A single-cell run
        }
        return lanes;
    };

    long long lapsChecked = 0, fillsChecked = 0, mismatches = 0;
    std::vector<SweepLane> lanes = randomLanes(60);
    SectorLapTable table(lanes);
    for (int round = 0; round < 3; ++round) {
        const int n = static_cast<int>(lanes.size());
        for (int q = 0; q < 4000; ++q) {
            const int first = static_cast<int>(rng() % n);
            const int bound = first + static_cast<int>(rng() % (n - first));
            // A limit on a real lap from this first lane, one float step either side of it, or random
            const float onLap = estimateSectorLength(lanes, first, first + static_cast<int>(rng() % (n - first)));
            float limit = onLap;
            switch (rng() % 4) {
            case 0: break;
            case 1: limit = std::nextafter(onLap, 0.0f); break;
            case 2: limit = std::nextafter(onLap, std::numeric_limits<float>::max()); break;
            default: limit = std::uniform_real_distribution<float>(0.0f, 20000.0f)(rng); break;
            }
            if (table.greedyLast(first, bound, limit) != greedyRef(lanes, first, bound, limit)) ++mismatches;
            ++fillsChecked;
            if (table.lap(first, bound) != estimateSectorLength(lanes, first, bound)) ++mismatches;
            ++lapsChecked;
        }
        // Infinite and NaN limits take the plain path and must still agree
        for (float limit : { std::numeric_limits<float>::infinity(), std::numeric_limits<float>::quiet_NaN() }) {
            for (int first = 0; first < n; first += 7) {
                if (table.greedyLast(first, n - 1, limit) != greedyRef(lanes, first, n - 1, limit)) ++mismatches;
                ++fillsChecked;
            }
        }
        // The lanes change (same count, as after joining loops): reset and carry on with the new ones
        lanes = randomLanes(n);
        table.reset();
    }
    std::printf("  (%lld laps and %lld greedy fills compared; %lld differ)\n", lapsChecked, fillsChecked, mismatches);
    CHECK(mismatches == 0);
}

// The planner's sector split (now read from a SectorLapTable) must give exactly the sectors the original
// split gives on the plan's final lanes: for each lane group, with the number of sectors the plan gave it, the
// original min-max split (binary search on the lap, greedy fill lane by lane, kept here unchanged) must cut the
// group at the same lanes. Run on the exactness fires and on twelve separate spot fires of different sizes: more
// loops than MAX_SWEEP_DRONES, so the planner joins small loops (moving lanes) and resets its lap table. (After a
// join every group has exactly one drone, so today the split there is the whole group either way; the case is
// covered so that a future change to the sizing rules can't use stale laps unnoticed.)
TEST(sweep_sectors_match_the_original_split) {
    using Range = std::pair<int, int>;
    // Original splitGroup (before SectorLapTable), word for word apart from the name
    auto splitRef = [](const std::vector<SweepLane>& lanes, int first, int last, int parts) {
        const int n = last - first + 1;
        parts = std::max(1, std::min(parts, n));
        auto fill = [&](float T) {
            std::vector<Range> out;
            int s = first;
            while (s <= last) {
                SectorLapBuilder b(lanes, s);
                while (b.last() + 1 <= last && b.lapWithNext() <= T) b.addNext();
                out.emplace_back(s, b.last());
                s = b.last() + 1;
            }
            return out;
        };
        float lo = 0.0f;
        for (int i = first; i <= last; ++i) lo = std::max(lo, estimateSectorLength(lanes, i, i));
        float hi = std::max(lo, estimateSectorLength(lanes, first, last));
        for (int it = 0; it < 40 && hi - lo > 0.5f; ++it) {
            const float mid = 0.5f * (lo + hi);
            if (static_cast<int>(fill(mid).size()) <= parts) hi = mid; else lo = mid;
        }
        std::vector<Range> out = fill(hi);
        while (static_cast<int>(out.size()) < parts) {
            int pick = -1;
            float worst = -1.0f;
            for (int i = 0; i < static_cast<int>(out.size()); ++i) {
                if (out[i].second == out[i].first) continue;
                const float l = estimateSectorLength(lanes, out[i].first, out[i].second);
                if (l > worst) { worst = l; pick = i; }
            }
            if (pick < 0) break;
            const Range r = out[pick];
            const int mid = (r.first + r.second) / 2;
            out[pick] = Range(r.first, mid);
            out.insert(out.begin() + pick + 1, Range(mid + 1, r.second));
        }
        return out;
    };

    std::vector<std::vector<std::uint8_t>> fires = exactnessFires();
    // 4 x 3 spot fires, 300 cells (3 km) apart: twelve perimeter loops, of different sizes (radius 2 to 13 cells)
    // so that no two loops have the same lanes
    std::vector<std::uint8_t> spots(N, 0);
    int radius = 2;
    for (int dx : { -450, -150, 150, 450 })
        for (int dy : { -300, 0, 300 }) {
            const auto d = discFire(CX + dx, CY + dy, radius++);
            for (int i = 0; i < N; ++i) spots[i] |= d[i];
        }
    fires.push_back(spots);

    int groupsChecked = 0, differ = 0;
    bool joined = false;
    for (const auto& fire : fires) {
        SweepSetup s(fire);
        const std::vector<SweepLane>& lanes = s.plan.lanes;
        if (lanes.empty()) continue;
        // Lane groups (contiguous runs sharing a loop id) and the plan's sectors in each (contiguous runs of
        // lanes sharing a sector)
        std::vector<Range> groups;
        for (int i = 0; i < static_cast<int>(lanes.size()); ++i) {
            if (groups.empty() || lanes[i].loop != lanes[groups.back().first].loop) groups.emplace_back(i, i);
            else groups.back().second = i;
        }
        if (s.plan.perimeterLoops.size() > groups.size()) joined = true;   // Loops were joined
        for (const Range& g : groups) {
            std::vector<Range> planned;
            for (int i = g.first; i <= g.second; ++i) {
                if (planned.empty() || lanes[i].sector != lanes[planned.back().first].sector) planned.emplace_back(i, i);
                else planned.back().second = i;
            }
            differ += splitRef(lanes, g.first, g.second, static_cast<int>(planned.size())) != planned;
            ++groupsChecked;
        }
    }
    std::printf("  (%d lane groups checked; %d differ; loops joined: %s)\n", groupsChecked, differ, joined ? "yes" : "no");
    CHECK(differ == 0);
    CHECK(joined);   // The spot fires really make the planner join loops
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

// Re-plans that reuse an old plan's whole-grid arrays (SweepFleet::recycle) must give exactly the plans that fresh
// arrays give. The sequence grows the fire, moves it (so the old plan's window is partly outside the new one and
// stale values there must be cleared) and shrinks it again.
TEST(recycled_plan_buffers_match_fresh_plans) {
    const std::vector<std::vector<std::uint8_t>> fires = {
        discFire(CX, CY, 25), discFire(CX, CY, 60), bayFire(), discFire(CX + 120, CY - 60, 40), discFire(CX, CY, 20) };
    MissionPlan prevFresh, prevRecycled, spare;
    int compared = 0, differences = 0;
    for (const auto& fire : fires) {
        DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
        live.compute(fire);
        Perimeter perim(fire);
        SafePlanner sp(live);
        auto fleetFor = [](const MissionPlan& prev) {
            SweepFleet fleet;
            for (int i = 0; i < SWEEP_DRONES; ++i) { fleet.ids.push_back(i); fleet.positions.push_back(Vec2(0.5f, GRID_HEIGHT - 0.5f)); }
            if (!prev.label.empty()) fleet.previous = &prev;
            return fleet;
        };
        // Fresh arrays
        std::mt19937 rngA(2024);
        MissionPlan fresh = planMission(live, perim.distance(), sp, rngA, fleetFor(prevFresh));
        // Recycled arrays (the plan before the previous one, as the controller does)
        std::mt19937 rngB(2024);
        SweepFleet fleet = fleetFor(prevRecycled);
        fleet.recycle = &spare;
        MissionPlan recycled = planMission(live, perim.distance(), sp, rngB, fleet);

        differences += fresh.label != recycled.label;
        differences += fresh.zone != recycled.zone;
        differences += fresh.owner != recycled.owner;
        differences += fresh.sectorDrone != recycled.sectorDrone;
        auto samePoints = [](const std::vector<Vec2>& a, const std::vector<Vec2>& b) {
            if (a.size() != b.size()) return false;
            for (size_t i = 0; i < a.size(); ++i) if (a[i].x != b[i].x || a[i].y != b[i].y) return false;
            return true;
        };
        differences += !samePoints(fresh.acoTour, recycled.acoTour);
        differences += fresh.sweepRoutes.size() != recycled.sweepRoutes.size();
        for (size_t r = 0; r < std::min(fresh.sweepRoutes.size(), recycled.sweepRoutes.size()); ++r)
            differences += !samePoints(fresh.sweepRoutes[r], recycled.sweepRoutes[r]);
        ++compared;

        prevFresh = std::move(fresh);
        spare = std::move(prevRecycled);       // The plan before last becomes the spare, as in MissionController
        prevRecycled = std::move(recycled);
    }
    std::printf("  (%d re-plans compared, %d differences)\n", compared, differences);
    CHECK(differences == 0);
    CHECK(spare.label.empty() || spare.label.size() == static_cast<size_t>(N));
}

// Only the ACO region (the outer zone, perimeter line out to SEARCH_DISTANCE) is split into Voronoi cells. The
// inner band, which the sweep lanes cover, gets no Voronoi label, and every point of interest is in the outer zone.
TEST(voronoi_covers_only_the_aco_region) {
    for (const auto& fire : { discFire(CX, CY, 30), bayFire() }) {
        SweepSetup s(fire);
        const MissionPlan& plan = s.plan;
        int inner = 0, outer = 0, wrong = 0;
        for (int i = 0; i < N; ++i) {
            inner += plan.zone[i] == 1;
            outer += plan.zone[i] == 2;
            // Labelled exactly where the cell is in the outer zone, with a valid region index
            const bool labelled = plan.label[i] >= 0;
            if (labelled != (plan.zone[i] == 2)) ++wrong;
            if (labelled && plan.label[i] >= static_cast<int>(plan.regions.size())) ++wrong;
        }
        int poisOutside = 0, regionsUsed = 0;
        for (const Region& rg : plan.regions) {
            if (rg.outerCount == 0) continue;
            ++regionsUsed;
            poisOutside += plan.zone[at(static_cast<int>(rg.poi.x), static_cast<int>(rg.poi.y))] != 2;
        }
        std::printf("  (%d inner-band cells unlabelled, %d outer cells in %d Voronoi regions)\n", inner, outer, regionsUsed);
        CHECK(inner > 0 && outer > 0);    // Both zones are there
        CHECK(wrong == 0);
        CHECK(poisOutside == 0);
        CHECK(regionsUsed >= MIN_SITES);
        CHECK(plan.acoTour.size() >= 3);
    }
}

// The ACO routing only changes when the perimeter does. While the fire grows inside the perimeter (no section
// has moved) a re-plan keeps the same regions, points of interest and tour; once a section moves, the tour is
// planned afresh. The kept tour must still be fire-safe.
TEST(aco_tour_changes_only_when_the_perimeter_moves) {
    // Identical tours, point for point
    auto sameTour = [](const std::vector<Vec2>& a, const std::vector<Vec2>& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) if (a[i].x != b[i].x || a[i].y != b[i].y) return false;
        return true;
    };
    // Runs frames until the controller re-plans (at most one re-plan interval); true if it did
    auto runToReplan = [](MissionController& m) {
        for (int f = 0; f <= REPLAN_INTERVAL; ++f) {
            m.update();
            if (m.takePlanChanged()) return true;
        }
        return false;
    };

    std::vector<std::uint8_t> fire = discFire(CX, CY, 25);
    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
    live.compute(fire);
    Perimeter perim(fire);
    MissionController mission(live, perim.distance(), Vec2(2.5f, GRID_HEIGHT - 2.5f));
    mission.markFireChanged();
    CHECK(runToReplan(mission));
    const std::vector<Vec2> first = mission.plan().acoTour;
    CHECK(!mission.plan().acoKept);   // First plan: nothing to keep
    CHECK(first.size() >= 3);

    // --- The fire grows 150 m, well short of the 600 m that moves a perimeter section ---
    fire = discFire(CX, CY, 40);
    live.compute(fire);
    CHECK(!perim.update(fire));       // The perimeter holds
    mission.markFireChanged();
    CHECK(runToReplan(mission));      // The sweep side still re-plans for the bigger fire...
    CHECK(mission.plan().acoKept);    // ...but the ACO side is carried over
    CHECK(sameTour(mission.plan().acoTour, first));

    // Carrying the tour over gives exactly what planning afresh would: the Voronoi only covers the outer zone,
    // which the fire growing inside the perimeter does not change
    {
        SafePlanner sp(live);
        std::mt19937 rng(2024);
        const MissionPlan fresh = planMission(live, perim.distance(), sp, rng);
        std::printf("  (tour of %zu waypoints kept; a fresh plan %s it)\n", first.size(),
                    sameTour(fresh.acoTour, first) ? "matches" : "differs from");
        CHECK(sameTour(fresh.acoTour, first));
    }

    // The kept tour is still fire-safe: every leg, including the wrap back to the start, clears the hard no-fly zone
    {
        SafePlanner sp(live);
        const std::vector<Vec2>& t = mission.plan().acoTour;
        int unsafe = 0;
        for (size_t i = 0; i < t.size(); ++i) unsafe += !sp.lineClear(t[i], t[(i + 1) % t.size()]);
        CHECK(unsafe == 0);
    }

    // --- The fire grows 700 m: past the trigger, so perimeter sections move and the tour is planned afresh ---
    fire = discFire(CX, CY, 95);
    live.compute(fire);
    CHECK(perim.update(fire));
    mission.markFireChanged();
    CHECK(runToReplan(mission));
    CHECK(!mission.plan().acoKept);
    CHECK(!sameTour(mission.plan().acoTour, first));
}

// While no perimeter section has moved, a re-plan copies the previous plan's perimeter line instead of tracing
// it again (it depends only on the committed-perimeter field, which is unchanged). The copy must be exactly
// what tracing gives, and the plan built on it (lanes, routes, ownership, tour) exactly the plan built on a fresh
// trace. If the field has changed, or the previous plan has no trace of it, the line is traced afresh.
TEST(perimeter_line_reused_only_while_the_perimeter_holds) {
    auto samePoints = [](const std::vector<Vec2>& a, const std::vector<Vec2>& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) if (a[i].x != b[i].x || a[i].y != b[i].y) return false;
        return true;
    };
    auto sameLoops = [&](const std::vector<std::vector<Vec2>>& a, const std::vector<std::vector<Vec2>>& b) {
        if (a.size() != b.size()) return false;
        for (size_t l = 0; l < a.size(); ++l) if (!samePoints(a[l], b[l])) return false;
        return true;
    };
    auto fleetAfter = [](const MissionPlan& prev, bool perimeterUnchanged) {
        SweepFleet fleet;
        for (int i = 0; i < SWEEP_DRONES; ++i) { fleet.ids.push_back(i); fleet.positions.push_back(Vec2(0.5f, GRID_HEIGHT - 0.5f)); }
        fleet.previous = &prev;
        fleet.perimeterUnchanged = perimeterUnchanged;
        return fleet;
    };

    // First plan, on a small fire
    std::vector<std::uint8_t> fire = discFire(CX, CY, 25);
    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
    live.compute(fire);
    Perimeter perim(fire);
    MissionPlan first;
    {
        SafePlanner sp(live);
        std::mt19937 rng(2024);
        first = planMission(live, perim.distance(), sp, rng);
    }
    CHECK(first.perimeterLoopsVersion == perim.distance().version());
    CHECK(!first.perimeterLoops.empty());

    // The fire grows inside the perimeter (no section moves): the line is reused...
    fire = discFire(CX, CY, 40);
    live.compute(fire);
    CHECK(!perim.update(fire));
    SafePlanner sp(live);
    std::mt19937 rngA(99);
    const MissionPlan reused = planMission(live, perim.distance(), sp, rngA, fleetAfter(first, true));
    const auto traced = tracePerimeterLine(perim.distance());
    CHECK(sameLoops(reused.perimeterLoops, traced));
    CHECK(sameLoops(reused.perimeterLoops, first.perimeterLoops));

    // ...and the plan is the one a fresh trace gives: same previous plan, but marked as having no trace
    MissionPlan untraced = first;
    untraced.perimeterLoopsVersion = 0;
    untraced.perimeterLoops.clear();
    std::mt19937 rngB(99);
    const MissionPlan fresh = planMission(live, perim.distance(), sp, rngB, fleetAfter(untraced, true));
    CHECK(sameLoops(fresh.perimeterLoops, traced));
    CHECK(sameLanes(reused.lanes, fresh.lanes));
    CHECK(reused.owner == fresh.owner);
    CHECK(reused.sectorDrone == fresh.sectorDrone);
    CHECK(reused.sweepRoutes.size() == fresh.sweepRoutes.size());
    for (size_t r = 0; r < std::min(reused.sweepRoutes.size(), fresh.sweepRoutes.size()); ++r)
        CHECK(samePoints(reused.sweepRoutes[r], fresh.sweepRoutes[r]));
    CHECK(samePoints(reused.acoTour, fresh.acoTour));

    // The perimeter moves: even if a caller wrongly says it has not, the version check traces the new line
    fire = discFire(CX, CY, 95);
    live.compute(fire);
    CHECK(perim.update(fire));
    SafePlanner sp2(live);
    std::mt19937 rngC(99);
    const MissionPlan moved = planMission(live, perim.distance(), sp2, rngC, fleetAfter(reused, true));
    CHECK(sameLoops(moved.perimeterLoops, tracePerimeterLine(perim.distance())));
    CHECK(!sameLoops(moved.perimeterLoops, first.perimeterLoops));
    CHECK(moved.perimeterLoopsVersion == perim.distance().version());
    std::printf("  (%zu-point line reused; plan identical to a fresh trace; moved perimeter traced afresh)\n",
                traced.empty() ? size_t(0) : traced[0].size());
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
