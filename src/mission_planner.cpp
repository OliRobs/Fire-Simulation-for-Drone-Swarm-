#include "mission_planner.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include "aco.h"
#include "boustrophedon.h"
#include "perimeter.h"
#include "voronoi.h"

// ============================================================================
// Mission Planner Engine
//
// Builds the complete multi-drone survey plan based on:
// 1. Live Fire Standoff (Dynamic No-Fly Safety Zone around active fire).
// 2. Committed Perimeter Zones:
//    - Inner band (between the standoff and the perimeter line) -> sweep lanes at right angles to the fire,
//      split into contiguous sectors along the perimeter line, one per sweep drone.
//    - Outer Zone (Between perimeter and search distance bound) -> ACO Perimeter Loop.
// 3. Voronoi Spatial Partitioning (k-means / Lloyd relaxation) to split the outer zone (the ACO region only;
//    the inner band is left to the sweep) into regions that give the ACO points of interest (POIs).
// ============================================================================

namespace {
    // ------------------------------------------------------------------------
    // Fire-Safe Closed Loop
    // ------------------------------------------------------------------------
    // Connects waypoints using A* pathfinding whenever straight-line flight intersects the active fire no-fly
    // zone. Waypoints that can only be reached by crossing the fire are dropped, so a route never contains an
    // unsafe leg. The last leg back to the start closes the loop.
    std::vector<Vec2> safeLoop(SafePlanner& planner, const std::vector<Vec2>& wp) {
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
    }

    // A contiguous run [first, last] of plan.lanes (inclusive)
    using LaneRange = std::pair<int, int>;

    // ------------------------------------------------------------------------
    // Split One Lane Group into Sectors
    // ------------------------------------------------------------------------
    // Splits lanes[first..last] into at most `parts` contiguous sectors so that the longest estimated lap is as
    // short as possible: binary-search the lap length, filling sectors greedily (each takes lanes while its
    // estimated lap still fits). A sector's estimate only grows as lanes are added (triangle inequality on the
    // return leg), so the greedy fill is exact for a given limit. Always returns exactly
    // min(parts, number of lanes) sectors: if the best split needs fewer, the largest ones are cut in two.
    std::vector<LaneRange> splitGroup(const std::vector<SweepLane>& lanes, int first, int last, int parts) {
        const int n = last - first + 1;
        parts = std::max(1, std::min(parts, n));

        // Greedy fill for limit T; returns the sectors (may be more than `parts` if T is too small)
        auto fill = [&](float T) {
            std::vector<LaneRange> out;
            int s = first;
            while (s <= last) {
                // Grow the sector one lane at a time while its lap still fits. The builder keeps the running
                // length, so each test costs one lane rather than re-adding the whole sector (same result
                // as estimateSectorLength(lanes, s, e + 1), bit for bit)
                SectorLapBuilder b(lanes, s);
                while (b.last() + 1 <= last && b.lapWithNext() <= T) b.addNext();
                out.emplace_back(s, b.last());
                s = b.last() + 1;
            }
            return out;
        };

        // The lap can't be shorter than the longest single-lane sector, nor longer than the whole group
        float lo = 0.0f;
        for (int i = first; i <= last; ++i) lo = std::max(lo, estimateSectorLength(lanes, i, i));
        float hi = std::max(lo, estimateSectorLength(lanes, first, last));
        for (int it = 0; it < 40 && hi - lo > 0.5f; ++it) {
            const float mid = 0.5f * (lo + hi);
            if (static_cast<int>(fill(mid).size()) <= parts) hi = mid; else lo = mid;
        }
        std::vector<LaneRange> out = fill(hi);

        // Too few sectors (spare drones): cut the longest multi-lane sector in two until there are `parts`
        while (static_cast<int>(out.size()) < parts) {
            int pick = -1;
            float worst = -1.0f;
            for (int i = 0; i < static_cast<int>(out.size()); ++i) {
                if (out[i].second == out[i].first) continue;
                const float l = estimateSectorLength(lanes, out[i].first, out[i].second);
                if (l > worst) { worst = l; pick = i; }
            }
            if (pick < 0) break;
            const LaneRange r = out[pick];
            const int mid = (r.first + r.second) / 2;
            out[pick] = LaneRange(r.first, mid);
            out.insert(out.begin() + pick + 1, LaneRange(mid + 1, r.second));
        }
        return out;
    }

    // Longest estimated lap among some sectors (cells)
    float longestEstimate(const std::vector<SweepLane>& lanes, const std::vector<LaneRange>& sectors) {
        float worst = 0.0f;
        for (const LaneRange& r : sectors) worst = std::max(worst, estimateSectorLength(lanes, r.first, r.second));
        return worst;
    }

    // ------------------------------------------------------------------------
    // Lane Groups
    // ------------------------------------------------------------------------
    // The contiguous ranges of plan.lanes that share a loop id (lanes are kept grouped by loop)
    std::vector<LaneRange> laneGroups(const std::vector<SweepLane>& lanes) {
        std::vector<LaneRange> groups;
        for (int i = 0; i < static_cast<int>(lanes.size()); ++i) {
            if (groups.empty() || lanes[i].loop != lanes[groups.back().first].loop) groups.emplace_back(i, i);
            else groups.back().second = i;
        }
        return groups;
    }

    // ------------------------------------------------------------------------
    // Join Small Loops onto Others
    // ------------------------------------------------------------------------
    // A sector never spans two loops, so with more loops than sweep drones (several spot fires, say) the smallest
    // loops are joined to their nearest neighbour: the small loop's lanes are moved, as one block, to just after
    // the lane of another loop nearest to it. The transit between them then counts towards that sector's lap
    // (as the step between two lanes). Repeats until there are at most maxGroups groups.
    void joinSmallLoops(std::vector<SweepLane>& lanes, int maxGroups) {
        maxGroups = std::max(1, maxGroups);
        for (;;) {
            const std::vector<LaneRange> groups = laneGroups(lanes);
            if (static_cast<int>(groups.size()) <= maxGroups) return;

            // Smallest loop by estimated length of flying it all
            int small = 0;
            float smallLen = std::numeric_limits<float>::max();
            for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
                const float l = estimateSectorLength(lanes, groups[g].first, groups[g].second);
                if (l < smallLen) { smallLen = l; small = g; }
            }
            const LaneRange sr = groups[small];

            // Nearest lane in any other loop (seed to seed)
            int host = -1;
            float hd = std::numeric_limits<float>::max();
            for (int i = 0; i < static_cast<int>(lanes.size()); ++i) {
                if (i >= sr.first && i <= sr.second) continue;
                for (int j = sr.first; j <= sr.second; ++j) {
                    const float d = sqDist(lanes[i].seed, lanes[j].seed);
                    if (d < hd) { hd = d; host = i; }
                }
            }
            const int hostLoop = lanes[host].loop;

            // Move the block after the host lane (and after any gap-fill lanes attached to it)
            std::vector<SweepLane> block(lanes.begin() + sr.first, lanes.begin() + sr.second + 1);
            for (SweepLane& l : block) l.loop = hostLoop;
            lanes.erase(lanes.begin() + sr.first, lanes.begin() + sr.second + 1);
            int at = host < sr.first ? host + 1 : host - (sr.second - sr.first + 1) + 1;
            while (at < static_cast<int>(lanes.size()) && lanes[at].gapFill) ++at;
            lanes.insert(lanes.begin() + at, block.begin(), block.end());
        }
    }

    // ------------------------------------------------------------------------
    // Share Drones Between Lane Groups
    // ------------------------------------------------------------------------
    // Gives each group one drone, then each further drone to the group whose longest sector is currently the
    // longest. For a min-max goal with laps that only shrink as drones are added, this greedy is optimal.
    // A group never gets more drones than it has lanes.
    std::vector<int> shareDrones(const std::vector<SweepLane>& lanes, const std::vector<LaneRange>& groups, int k) {
        const int g = static_cast<int>(groups.size());
        std::vector<int> share(g, 1);
        std::vector<float> worst(g);
        for (int i = 0; i < g; ++i) worst[i] = longestEstimate(lanes, splitGroup(lanes, groups[i].first, groups[i].second, 1));
        for (int left = k - g; left > 0; --left) {
            int pick = -1;
            for (int i = 0; i < g; ++i) {
                if (share[i] >= groups[i].second - groups[i].first + 1) continue;   // One lane each already
                if (pick < 0 || worst[i] > worst[pick]) pick = i;
            }
            if (pick < 0) break;
            ++share[pick];
            worst[pick] = longestEstimate(lanes, splitGroup(lanes, groups[pick].first, groups[pick].second, share[pick]));
        }
        return share;
    }

    // The sectors for a given share of drones per group, in sequence order
    std::vector<LaneRange> sectorsFor(const std::vector<SweepLane>& lanes, const std::vector<LaneRange>& groups,
                                      const std::vector<int>& share) {
        std::vector<LaneRange> out;
        for (size_t i = 0; i < groups.size(); ++i) {
            const std::vector<LaneRange> s = splitGroup(lanes, groups[i].first, groups[i].second, share[i]);
            out.insert(out.end(), s.begin(), s.end());
        }
        return out;
    }

    // ------------------------------------------------------------------------
    // Drones Needed for One Lane Group
    // ------------------------------------------------------------------------
    // Fewest sectors this group can be split into with every estimated lap within `limit` cells (capped at
    // `most`, or the number of lanes if smaller)
    int dronesNeeded(const std::vector<SweepLane>& lanes, LaneRange group, float limit, int most) {
        const int n = group.second - group.first + 1;
        most = std::max(1, std::min(most, n));
        for (int k = 1; k < most; ++k) {
            if (longestEstimate(lanes, splitGroup(lanes, group.first, group.second, k)) <= limit) return k;
        }
        return most;
    }

    // ------------------------------------------------------------------------
    // Match Sectors to Drones (stable ownership)
    // ------------------------------------------------------------------------
    // Returns, for each sector, the id of the drone that should fly it, or -1 if no drone is left for it (a new
    // drone must be launched). First by shared land: count the band cells each new sector shares with each
    // drone's land in the previous plan, then repeatedly give the pair with the most shared cells to each other
    // (greedy). Sectors are cut along the perimeter line in a fixed direction, so a drone that keeps its side of
    // the fire keeps most of its sector. Sectors and drones still unmatched after that (first plan, new land) are
    // paired by distance: each drone to the nearest unclaimed sector start, nearest pairs first.
    std::vector<int> matchSectors(const MissionPlan& plan, const std::vector<std::vector<Vec2>>& routes,
                                  const SweepFleet& fleet, const GridRect& area) {
        const int S = static_cast<int>(routes.size()), D = static_cast<int>(fleet.ids.size());
        std::vector<int> result(S, -1);
        std::vector<char> droneUsed(D, 0);

        struct Pair { double score; int sector, drone; };
        std::vector<Pair> pairs;

        // --- By shared cells with the previous plan ---
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

        // --- The rest by distance from each drone to the nearest waypoint of each sector ---
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

        // --- Sectors without a route (no lanes at all): any drone still free ---
        for (int s = 0; s < S; ++s) {
            if (result[s] >= 0) continue;
            for (int d = 0; d < D; ++d) if (!droneUsed[d]) { result[s] = fleet.ids[d]; droneUsed[d] = 1; break; }
        }
        return result;
    }
}

// ----------------------------------------------------------------------------
// Plan for the starting fleet
// ----------------------------------------------------------------------------
MissionPlan planMission(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt,
    SafePlanner& planner, std::mt19937& rng) {
    SweepFleet fleet;
    for (int i = 0; i < SWEEP_DRONES; ++i) {
        fleet.ids.push_back(i);
        fleet.positions.push_back(Vec2(0.5f, GRID_HEIGHT - 0.5f));
    }
    return planMission(liveEdt, perimEdt, planner, rng, fleet);
}

MissionPlan planMission(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt,
    SafePlanner& planner, std::mt19937& rng, const SweepFleet& fleet) {
    constexpr int N = GRID_WIDTH * GRID_HEIGHT;
    constexpr double innerHi = (PERIMETER_DISTANCE + 0.5) * (PERIMETER_DISTANCE + 0.5);
    constexpr double outerHi = (SEARCH_DISTANCE + 0.5) * (SEARCH_DISTANCE + 0.5);
    const int numSweep = static_cast<int>(fleet.ids.size());

    MissionPlan plan;
    plan.label.assign(N, -1);
    plan.zone.assign(N, 0);   // 1 = inner zone, 2 = outer zone, 0 = unmapped / no-fly
    plan.sweepRoutes.assign(numSweep, {});
    plan.sectorDrone = fleet.ids;   // Until sectors are made (empty search area: each drone keeps an empty one)

    // ------------------------------------------------------------------------
    // Step 1: Classify Grid Cells into Inner, Outer, or Invalid Zones
    // ------------------------------------------------------------------------
    std::vector<std::uint8_t>& zone = plan.zone;
    bool anyCells = false;     // Any surveyable cell at all (inner band or outer zone)
    std::vector<int> outer;    // Outer-zone cells (the ACO region), in scan order: the only cells the Voronoi splits

    // Only the perimeter transform's window can hold search-area cells: outside it every cell reads INF (for
    // the Perimeter's transform, more than SEARCH_DISTANCE + 1 from the committed fire), so beyond outerHi.
    // Scanning just the window, row by row, finds the same cells in the same order as scanning the whole grid.
    const GridRect win = perimEdt.window();
    plan.bounds = win;
    for (int y = win.y0; y < win.y1; ++y) {
        for (int i = y * GRID_WIDTH + win.x0, end = y * GRID_WIDTH + win.x1; i < end; ++i) {
            // Exclude cells inside the no-fly standoff zone around the active fire
            if (liveEdt.squaredDistance(i) < NO_FLY_SQ) continue;

            // Categorize valid cells based on squared distance from committed perimeter
            const double d2 = perimEdt.squaredDistance(i);
            if (d2 < innerHi) {
                zone[i] = 1; // Inner zone cell (swept by lanes, not split into Voronoi regions)
            }
            else if (d2 < outerHi) {
                zone[i] = 2; // Outer zone cell (ACO region)
                outer.push_back(i);
            }
            else {
                continue;    // Beyond max search boundary
            }
            anyCells = true;
        }
    }

    // Return empty plan if no surveyable cells exist
    if (!anyCells) return plan;

    // ------------------------------------------------------------------------
    // Step 2a: Keep the ACO side of the previous plan while the perimeter holds
    // ------------------------------------------------------------------------
    // The Voronoi regions only supply the ACO points of interest, and cover only the outer zone (perimeter
    // line out to SEARCH_DISTANCE). If the committed perimeter has not moved, the outer zone is the same land
    // (the live fire is still at least PERIMETER_UPDATE_DISTANCE inside the perimeter line, well clear of it),
    // so the regions, their points of interest and the tour order are carried over unchanged and the ACO
    // routing stays put. Rebuilding them would give the same result (same cells, same seed); carrying them
    // over skips the Voronoi and the ant colony, and guarantees the tour holds. Every outer cell must have been
    // labelled in the previous plan; if not (e.g. a caller set the flag wrongly) the regions are rebuilt.
    const MissionPlan* prev = fleet.previous;
    plan.acoKept = fleet.perimeterUnchanged && prev && !prev->regions.empty()
                   && prev->label.size() == static_cast<size_t>(N);
    if (plan.acoKept) {
        for (int idx : outer) {
            if (prev->label[idx] < 0) { plan.acoKept = false; break; }
        }
    }
    if (plan.acoKept) {
        for (int idx : outer) plan.label[idx] = prev->label[idx];
        plan.regions = prev->regions;   // Same outer-zone cells, so the same points of interest
    }
    else {
        // ------------------------------------------------------------------------
        // Step 2: Build Voronoi Spatial Regions (Lloyd's Relaxation)
        // ------------------------------------------------------------------------
        // Only the outer zone (the ACO region) is split into regions; the inner band is covered by the sweep
        // lanes and gets no label. Region count 'S' from the outer zone's area, clamped between MIN_SITES and
        // MAX_SITES (and never more regions than cells)
        const int S = std::min(std::clamp(static_cast<int>(outer.size()) / TARGET_REGION_CELLS, MIN_SITES, MAX_SITES),
            static_cast<int>(outer.size()));

        const Voronoi voronoi = buildVoronoi(outer, S, rng);

        // Label each outer-zone cell with its assigned Voronoi site index
        for (size_t i = 0; i < outer.size(); ++i) {
            plan.label[outer[i]] = voronoi.assign[i];
        }

        // ------------------------------------------------------------------------
        // Step 3: Outer-Zone Statistics per Region and their Points of Interest
        // ------------------------------------------------------------------------
        plan.regions.assign(S, Region());
        for (int idx : outer) {
            Region& rg = plan.regions[plan.label[idx]];
            ++rg.outerCount;
            rg.outerSum += cellCentre(idx);
        }
        for (Region& rg : plan.regions) {
            if (rg.outerCount > 0) rg.outerCentroid = rg.outerSum / static_cast<float>(rg.outerCount);
        }

        // Snap outer Point of Interest (POI) to the nearest actual cell in the outer region
        for (int idx : outer) {
            Region& rg = plan.regions[plan.label[idx]];
            const float dd = sqDist(cellCentre(idx), rg.outerCentroid);
            if (dd < rg.poiDist) {
                rg.poiDist = dd;
                rg.poi = cellCentre(idx);
            }
        }
    }   // End of Steps 2-3 (new regions)

    // ------------------------------------------------------------------------
    // Step 4: Sweep Lanes along the Perimeter Line
    // ------------------------------------------------------------------------
    // Trace the 1.2 km line into loops, seed a lane at most every 120 m along each, run each lane straight in
    // towards the fire, then add short gap-fill lanes where lanes spread apart in concave bays.
    plan.perimeterLoops = tracePerimeterLine(perimEdt);
    plan.lanes = layLanes(plan.perimeterLoops, perimEdt, zone);
    fillCoverageGaps(plan.lanes, perimEdt, zone, win);

    // ------------------------------------------------------------------------
    // Step 5: Size the Sweep Fleet and Split the Lanes into Sectors
    // ------------------------------------------------------------------------
    // Lap limits in cells of flying: the hard cap, and the 90 % target the straight-line estimates are sized to
    // (the other 10 % allows for fire detours, which only show up once the routes are built)
    const float cellsPerSecond = DRONE_SPEED_MPS / SIM_CELL_SIZE_M;
    const float capCells = MAX_SWEEP_LAP_S * cellsPerSecond;
    const float targetCells = capCells * SWEEP_LAP_MARGIN;

    std::vector<LaneRange> sectors;
    if (!plan.lanes.empty()) {
        // 1-2. Estimate: per loop, the fewest contiguous sectors whose estimated laps (lanes, steps between
        //      lanes, straight return leg) all fit in the target; their sum is the drones needed
        std::vector<LaneRange> groups = laneGroups(plan.lanes);
        int needed = 0;
        for (const LaneRange& g : groups) needed += dronesNeeded(plan.lanes, g, targetCells, MAX_SWEEP_DRONES);

        // Drones to plan for: never fewer than the fleet has now (shrinking is the controller's call, made
        // slowly) and within the configured range
        int k = std::clamp(std::max(needed, numSweep), MIN_SWEEP_DRONES, MAX_SWEEP_DRONES);

        // A sector never spans two loops; with more loops than drones, join the smallest onto their neighbours
        if (static_cast<int>(groups.size()) > k) {
            joinSmallLoops(plan.lanes, k);
            groups = laneGroups(plan.lanes);
        }

        // Builds the real fire-safe routes for one share of drones per loop, and times every lap
        struct Split {
            std::vector<LaneRange> sectors;
            std::vector<std::vector<Vec2>> routes;
            std::vector<float> laps;
            int worst = -1;           // Sector with the longest lap
        };
        auto build = [&](const std::vector<int>& share) {
            Split sp;
            sp.sectors = sectorsFor(plan.lanes, groups, share);
            for (int s = 0; s < static_cast<int>(sp.sectors.size()); ++s) {
                sp.routes.push_back(safeLoop(planner, sectorWaypoints(plan.lanes, sp.sectors[s].first, sp.sectors[s].second)));
                sp.laps.push_back(flightSeconds(closedRouteLength(sp.routes.back())));
                if (sp.worst < 0 || sp.laps[s] > sp.laps[sp.worst]) sp.worst = s;
            }
            return sp;
        };
        // Which loop a sector belongs to (groups and sectors are both in sequence order)
        auto groupOf = [&](const LaneRange& sector) {
            int g = 0;
            while (sector.first > groups[g].second) ++g;
            return g;
        };

        // 3. Split: share the k drones between loops (worst loop first) and split each loop min-max
        std::vector<int> share = shareDrones(plan.lanes, groups, k);

        // 4. Check: build the real fire-safe routes and time them. If a lap is over the cap, add a drone to that
        //    lap's loop and split again (at most 3 extra passes)
        Split split = build(share);
        for (int pass = 0; pass < 3; ++pass) {
            const int total = std::accumulate(share.begin(), share.end(), 0);
            if (split.worst < 0 || split.laps[split.worst] <= MAX_SWEEP_LAP_S || total >= MAX_SWEEP_DRONES) break;
            const int g = groupOf(split.sectors[split.worst]);
            if (share[g] >= groups[g].second - groups[g].first + 1) break;   // Already one lane per drone
            ++share[g];
            needed = std::max(needed, total + 1);
            split = build(share);
        }

        // 5. Shrink check: only when the plan uses exactly the drones there are (no growth this time) and one
        //    could go without breaking MIN_SWEEP_DRONES or a loop's own drone. Take the drone from the loop where
        //    losing it lengthens the estimated laps least, then build and time that k - 1 split for real.
        const int total = std::accumulate(share.begin(), share.end(), 0);
        if (total == numSweep && total > MIN_SWEEP_DRONES) {
            int from = -1;
            float best = std::numeric_limits<float>::max();
            for (int g = 0; g < static_cast<int>(groups.size()); ++g) {
                if (share[g] <= 1) continue;
                const float l = longestEstimate(plan.lanes, splitGroup(plan.lanes, groups[g].first, groups[g].second, share[g] - 1));
                if (l < best) { best = l; from = g; }
            }
            if (from >= 0) {
                std::vector<int> fewer = share;
                --fewer[from];
                Split smaller = build(fewer);
                plan.sweepShrinkOk = smaller.worst >= 0 && smaller.laps[smaller.worst] <= SHRINK_MARGIN * MAX_SWEEP_LAP_S;
                // Third (SHRINK_REPLANS-th) re-plan in a row that passes: use the smaller split
                if (plan.sweepShrinkOk && fleet.shrinkStreak + 1 >= SHRINK_REPLANS) {
                    split = std::move(smaller);
                    plan.sweepShrank = true;
                }
            }
        }

        sectors = split.sectors;
        plan.sweepRoutes = split.routes;
        plan.sweepLapS = split.laps;
        plan.sweepNeeded = std::clamp(needed, MIN_SWEEP_DRONES, MAX_SWEEP_DRONES);
        for (float lap : plan.sweepLapS) plan.sweepOverLimit = plan.sweepOverLimit || lap > MAX_SWEEP_LAP_S;
    }
    else {
        // No lanes (the fire covers the whole band): one empty sector per drone, so each keeps its place
        sectors.assign(std::max(numSweep, 0), LaneRange(0, -1));
        plan.sweepRoutes.assign(sectors.size(), {});
        plan.sweepLapS.assign(sectors.size(), 0.0f);
        plan.sweepNeeded = std::clamp(numSweep, MIN_SWEEP_DRONES, MAX_SWEEP_DRONES);
    }
    for (int s = 0; s < static_cast<int>(sectors.size()); ++s) {
        for (int i = sectors[s].first; i <= sectors[s].second; ++i) plan.lanes[i].sector = s;
    }

    // Ownership: each band cell to its nearest lane's sector, then sectors to drones by shared land
    plan.owner = buildOwnerMap(plan.lanes, zone, win);
    plan.sectorDrone = matchSectors(plan, plan.sweepRoutes, fleet, win);

    // ------------------------------------------------------------------------
    // Step 6: Plan Outer Zone Ant Colony Optimization (ACO) Tour
    // ------------------------------------------------------------------------
    // Perimeter unchanged: the same points of interest in the same order as last time. Otherwise solve the TSP
    // over the outer POIs with the Elitist Ant System.
    std::vector<Vec2> tour;
    if (plan.acoKept) {
        tour = prev->acoOrder;
    }
    else {
        std::vector<Vec2> pois;
        for (const Region& rg : plan.regions) {
            if (rg.outerCount > 0) pois.push_back(rg.poi);
        }
        for (int i : solveACO(pois, rng)) tour.push_back(pois[i]);
    }
    plan.acoOrder = tour;

    // Make it fire-safe. This is redone on every re-plan, kept tour or not: the live fire may have grown across
    // a leg since. Where every leg is still clear, the same order gives exactly the same route as before.
    plan.acoTour = safeLoop(planner, tour);

    // Always fly the tour anticlockwise on screen (negative signed area in y-down grid coordinates, the same
    // sense as the perimeter line). The ant colony returns the cycle in either direction, so without this a
    // re-plan could reverse the tour and send every ACO drone back the way it came. Reversing the waypoint
    // order keeps the same fire-safe legs, flown the other way.
    if (plan.acoTour.size() >= 3 && loopSignedArea(plan.acoTour) > 0.0) {
        std::reverse(plan.acoTour.begin() + 1, plan.acoTour.end());   // Keep the same start point
    }

    // ACO drone count: one drone takes acoLoopS to fly the safe tour; n evenly spaced drones pass each point
    // every acoLoopS / n, so n = ceil(acoLoopS / ACO_REVISIT_S) keeps that within the revisit limit
    plan.acoLoopS = flightSeconds(closedRouteLength(plan.acoTour));
    const int acoWanted = static_cast<int>(std::ceil(plan.acoLoopS / ACO_REVISIT_S - 1e-6f));
    plan.acoNeeded = std::clamp(acoWanted, MIN_ACO_DRONES, MAX_ACO_DRONES);
    plan.acoOverLimit = acoWanted > MAX_ACO_DRONES;

    return plan;
}
