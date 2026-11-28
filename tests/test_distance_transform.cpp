// ============================================================================
// Tests: Distance transform (distance_transform.h): full, range-limited and live incremental fields
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <vector>

#include "distance_transform.h"
#include "fire_model.h"
#include "test_framework.h"
#include "test_support.h"

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
