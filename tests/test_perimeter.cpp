// ============================================================================
// Tests: Perimeter (perimeter.h): section updates and the ordered perimeter line
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "distance_transform.h"
#include "fire_model.h"
#include "viewport.h"
#include "perimeter.h"
#include "test_framework.h"
#include "test_support.h"

// Perimeter::update with the fire's rectangle (main.cpp) and without it (scanning) must move the perimeter
// identically: same ring, same distances everywhere the mission planner reads them
SLOW_TEST(perimeter_update_with_bounds_matches_scan) {
    std::vector<Cell> grid(N);
    FireSimulation sim(G, grid);
    sim.ignite(at(400, 400));
    for (int t = 0; t < 200; ++t) sim.step();
    Perimeter withBounds(G, sim.fireMask()), scanned(G, sim.fireMask());
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
        FireSimulation sim(G, grid);
        sim.ignite(at(CX, CY));
        for (int t = 0; t < 100; ++t) sim.step();
        Perimeter fast(G, sim.fireMask());
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
        Perimeter fast(G, mask);
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
// Perimeter line (ordered loops traced from the committed-fire distance field)
// ============================================================================
TEST(perimeter_line_round_fire_is_one_loop) {
    const int r = 30;
    Perimeter perim(G, discFire(400, 400, r));
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
    Perimeter perim(G, fire);
    const auto loops = tracePerimeterLine(perim.distance());
    CHECK(loops.size() == 2);
    for (const auto& l : loops) CHECK(loopSignedArea(l) < 0.0);
}
