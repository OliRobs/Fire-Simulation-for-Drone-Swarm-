// ============================================================================
// Tests: Planning algorithms: ACO tour (aco.h), Voronoi regions (voronoi.h), planned routes
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include "aco.h"
#include "distance_transform.h"
#include "mission_planner.h"
#include "perimeter.h"
#include "safe_planner.h"
#include "voronoi.h"
#include "test_framework.h"
#include "test_support.h"

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
    const Voronoi v = buildVoronoi(G, cells, 9, rng);
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
        const Voronoi v = buildVoronoi(G, cells, sites, rng);
        CHECK(static_cast<int>(v.sites.size()) == sites);
        int wrong = 0;
        for (size_t i = 0; i < cells.size(); ++i) {
            const Vec2 p = G.centre(cells[i]);
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
        for (size_t i = 0; i < n; ++i) pts[i] = G.centre(cells[i]);
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
            const Voronoi fast = buildVoronoi(G, *cells, sites, rngA);
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
    Perimeter perim(G, fire);
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
