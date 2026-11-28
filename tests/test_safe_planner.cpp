// ============================================================================
// Tests: Safe planner (safe_planner.h): fire-safe legs, A* and line-of-sight checks
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <initializer_list>
#include <limits>
#include <queue>
#include <random>
#include <vector>

#include "distance_transform.h"
#include "safe_planner.h"
#include "test_framework.h"
#include "test_support.h"

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
SLOW_TEST(astar_matches_original_planner) {
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
SLOW_TEST(lineclear_skip_matches_every_sample) {
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
