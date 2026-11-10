#include "mission_planner.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include "aco.h"
#include "boustrophedon.h"
#include "voronoi.h"

// ============================================================================
// Mission Planner Engine
//
// Builds the complete multi-drone survey plan based on:
// 1. Live Fire Standoff (Dynamic No-Fly Safety Zone around active fire).
// 2. Committed Perimeter Zones:
//    - Inner Zone (Between standoff and perimeter) -> Boustrophedon Lawn-mower Sweep.
//    - Outer Zone (Between perimeter and search distance bound) -> ACO Perimeter Loop.
// 3. Voronoi Spatial Partitioning (k-means / Lloyd relaxation) to split the domain
//    into sub-regions and extract Points of Interest (POIs).
// ============================================================================
MissionPlan planMission(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt,
    SafePlanner& planner, std::mt19937& rng) {
    constexpr int N = GRID_WIDTH * GRID_HEIGHT;
    constexpr double innerHi = (PERIMETER_DISTANCE + 0.5) * (PERIMETER_DISTANCE + 0.5);
    constexpr double outerHi = (SEARCH_DISTANCE + 0.5) * (SEARCH_DISTANCE + 0.5);

    MissionPlan plan;
    plan.label.assign(N, -1);
    plan.sweepRoutes.assign(SWEEP_DRONES, {});

    // ------------------------------------------------------------------------
    // Step 1: Classify Grid Cells into Inner, Outer, or Invalid Zones
    // ------------------------------------------------------------------------
    std::vector<std::uint8_t> zone(N, 0); // 1 = inner zone, 2 = outer zone, 0 = unmapped / no-fly
    std::vector<int> cells;

    for (int i = 0; i < N; ++i) {
        // Exclude cells inside the no-fly standoff zone around the active fire
        if (liveEdt.squaredDistance(i) < NO_FLY_SQ) continue;

        // Categorize valid cells based on squared distance from committed perimeter
        const double d2 = perimEdt.squaredDistance(i);
        if (d2 < innerHi) {
            zone[i] = 1; // Inner zone cell
        }
        else if (d2 < outerHi) {
            zone[i] = 2; // Outer zone cell
        }
        else {
            continue;    // Beyond max search boundary
        }
        cells.push_back(i);
    }

    // Return empty plan if no surveyable cells exist
    if (cells.empty()) return plan;

    // ------------------------------------------------------------------------
    // Step 2: Build Voronoi Spatial Regions (Lloyd's Relaxation)
    // ------------------------------------------------------------------------
    // Calculate optimal region count 'S' clamped between MIN_SITES and MAX_SITES
    const int S = std::min(std::clamp(static_cast<int>(cells.size()) / TARGET_REGION_CELLS, MIN_SITES, MAX_SITES),
        static_cast<int>(cells.size()));

    const Voronoi voronoi = buildVoronoi(cells, S, rng);

    // Label each valid cell with its assigned Voronoi site index
    for (size_t i = 0; i < cells.size(); ++i) {
        plan.label[cells[i]] = voronoi.assign[i];
    }

    // ------------------------------------------------------------------------
    // Step 3: Compute Per-Region Statistics and Extract Outer POIs
    // ------------------------------------------------------------------------
    plan.regions.assign(S, Region());
    Vec2 innerSumAll(0.0f, 0.0f);
    int innerAll = 0;

    // Accumulate cell coordinates per region split by zone type
    for (int idx : cells) {
        Region& rg = plan.regions[plan.label[idx]];
        const int x = idx % GRID_WIDTH, y = idx / GRID_WIDTH;
        const Vec2 p = cellCentre(idx);

        if (zone[idx] == 1) { // Inner zone stats & bounding boxes
            ++rg.innerCount;
            rg.innerSum += p;
            rg.minX = std::min(rg.minX, x); rg.maxX = std::max(rg.maxX, x);
            rg.minY = std::min(rg.minY, y); rg.maxY = std::max(rg.maxY, y);
            innerSumAll += p;
            ++innerAll;
        }
        else {              // Outer zone stats
            ++rg.outerCount;
            rg.outerSum += p;
        }
    }

    // Calculate centroids for inner and outer region components
    for (Region& rg : plan.regions) {
        if (rg.innerCount > 0) rg.innerCentroid = rg.innerSum / static_cast<float>(rg.innerCount);
        if (rg.outerCount > 0) rg.outerCentroid = rg.outerSum / static_cast<float>(rg.outerCount);
    }

    // Snap outer Point of Interest (POI) to the nearest actual cell in the outer region
    for (int idx : cells) {
        if (zone[idx] != 2) continue;
        Region& rg = plan.regions[plan.label[idx]];
        const float dd = sqDist(cellCentre(idx), rg.outerCentroid);
        if (dd < rg.poiDist) {
            rg.poiDist = dd;
            rg.poi = cellCentre(idx);
        }
    }

    // ------------------------------------------------------------------------
    // Step 4: Plan Inner Zone Boustrophedon Sweep Routes
    // ------------------------------------------------------------------------
    // Order regions angularly around the global fire center and partition them
    // into balanced workload sectors assigned across available sweep drones.
    std::vector<int> innerRegions;
    int totalInner = 0;
    for (int r = 0; r < S; ++r) {
        if (plan.regions[r].innerCount > 0) {
            innerRegions.push_back(r);
            totalInner += plan.regions[r].innerCount;
        }
    }

    const Vec2 fireCentre = innerAll > 0 ? innerSumAll / static_cast<float>(innerAll) : Vec2(0.0f, 0.0f);

    auto angleOf = [&](int r) {
        const Vec2 c = plan.regions[r].innerCentroid;
        return std::atan2(c.y - fireCentre.y, c.x - fireCentre.x);
        };

    // Sort inner regions by polar angle relative to fire center
    std::sort(innerRegions.begin(), innerRegions.end(), [&](int a, int b) {
        return angleOf(a) < angleOf(b);
        });

    // Partition sorted sectors and append boustrophedon sweep paths per drone
    int cumulative = 0;
    for (int r : innerRegions) {
        const Region& rg = plan.regions[r];
        const double mid = cumulative + rg.innerCount * 0.5;
        const int chunk = std::min(SWEEP_DRONES - 1, static_cast<int>(mid * SWEEP_DRONES / totalInner));
        cumulative += rg.innerCount;
        appendBoustrophedon(r, rg, plan.label, zone, plan.sweepRoutes[chunk]);
    }

    // ------------------------------------------------------------------------
    // Step 5: Plan Outer Zone Ant Colony Optimization (ACO) Tour
    // ------------------------------------------------------------------------
    std::vector<Vec2> pois;
    for (const Region& rg : plan.regions) {
        if (rg.outerCount > 0) pois.push_back(rg.poi);
    }

    // Solve TSP over outer POIs using Elitist Ant System
    std::vector<Vec2> tour;
    for (int i : solveACO(pois, rng)) tour.push_back(pois[i]);

    // ------------------------------------------------------------------------
    // Step 6: Convert All Flight Paths into Fire-Safe Routes
    // ------------------------------------------------------------------------
    // Connects waypoints using A* pathfinding whenever straight-line flight
    // intersects the active fire no-fly zone. Waypoints that can only be reached by
    // crossing the fire are dropped, so a route never contains an unsafe leg.
    auto safeLoop = [&](const std::vector<Vec2>& wp) -> std::vector<Vec2> {
        if (wp.size() < 2) return wp;
        std::vector<Vec2> out{ wp[0] };

        for (size_t i = 1; i < wp.size(); ++i) {
            // Plan from where the route actually got to (out.back()), which can differ from wp[i - 1]
            // if that waypoint was dropped or the leg stopped short of it
            const std::vector<Vec2> l = planner.leg(out.back(), wp[i]);
            if (l.empty()) continue;   // Cut off by fire: skip this waypoint
            out.insert(out.end(), l.begin(), l.end());
        }
        if (out.size() < 2) return out;

        // Connect route end back to start point to complete the closed loop. The final point of the
        // closing leg is out.front() itself, so it is left off (the drone wraps round to index 0).
        const std::vector<Vec2> closing = planner.leg(out.back(), out.front());
        if (!closing.empty()) {
            out.insert(out.end(), closing.begin(), closing.end() - 1);
        }
        else {
            // No safe way back to the start: fly back along the same route instead (out-and-back), so the
            // wrap from the last waypoint to the first is always a leg that was already checked as safe
            for (size_t i = out.size() - 1; i-- > 1;) out.push_back(out[i]);
        }
        return out;
        };

    for (auto& route : plan.sweepRoutes) route = safeLoop(route);
    plan.acoTour = safeLoop(tour);

    return plan;
}