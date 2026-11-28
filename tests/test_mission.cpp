// ============================================================================
// Tests: Drones and the mission controller (drone.h, mission_controller.h): end-to-end mission tests
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <initializer_list>
#include <limits>
#include <random>
#include <set>
#include <vector>

#include "distance_transform.h"
#include "drone.h"
#include "fire_model.h"
#include "mission_controller.h"
#include "mission_planner.h"
#include "perimeter.h"
#include "safe_planner.h"
#include "voronoi.h"
#include "test_framework.h"
#include "test_support.h"

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
SLOW_TEST(mission_drones_never_fly_over_fire) {
    std::vector<Cell> grid(N);
    FireSimulation sim(G, grid);                       // Default weather: wind from the SW
    // Several fires so drones have to work round them
    sim.ignite(at(250, 250)); sim.ignite(at(180, 300)); sim.ignite(at(320, 200));

    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT, LIVE_EDT_RANGE);   // Range-limited, as in main.cpp
    for (int t = 0; t < 400; ++t) sim.step();      // Let the fire get going before drones launch
    live.compute(sim.fireMask());
    Perimeter perim(G, sim.fireMask());
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
    Perimeter perim(G, fire);
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
    Perimeter perim(G, fire);
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
    Perimeter perim(G, fire);
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
    Perimeter perim(G, fire);
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
        Perimeter perim(G, fire);
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
    Perimeter perim(G, fire);
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
    Perimeter perim(G, fire);
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
