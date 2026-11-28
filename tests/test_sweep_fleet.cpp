// ============================================================================
// Tests: Sweep fleet sizing: sector laps, lap caps and stable sector ownership
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <limits>
#include <random>
#include <vector>

#include "boustrophedon.h"
#include "distance_transform.h"
#include "mission_planner.h"
#include "perimeter.h"
#include "safe_planner.h"
#include "test_framework.h"
#include "test_support.h"

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
