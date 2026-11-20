#pragma once
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include "types.h"

// One Voronoi region of the outer zone (the ACO region, perimeter line out to SEARCH_DISTANCE). Regions only
// supply the ACO points of interest, one per region; the inner band is not split into regions at all (the sweep
// is planned from lanes along the perimeter line instead, SweepLane).
struct Region {
    int outerCount = 0;                   // Cells of this region outside the perimeter line (outer zone)
    Vec2 outerSum{ 0.0f, 0.0f };
    Vec2 outerCentroid{ 0.0f, 0.0f };
    Vec2 poi{ 0.0f, 0.0f };               // Outer-zone point of interest (visited by the ACO tour)
    float poiDist = std::numeric_limits<float>::max();
};

// One sweep lane: a straight line from the perimeter line in towards the fire, at right angles to the nearest
// fire edge, clipped to the inner band. If the line crosses a gap in the band it has several runs.
struct SweepLane {
    std::vector<Vec2> pts;    // Run end points in inward order: run 1 start, run 1 end, run 2 start, ...
                              // (a run of a single cell has the same start and end point)
    Vec2 seed{ 0.0f, 0.0f };  // Where it was seeded on the perimeter line (gap-fill lanes: their outer end)
    int loop = 0;             // Lane group it belongs to: the perimeter-line loop it was seeded on
    int sector = -1;          // Sector (index into MissionPlan::sweepRoutes) that flies it
    bool gapFill = false;     // Added by gap fill rather than seeded on the perimeter line

    // Length flown along the lane, end to end (including any gap between its runs), in cells
    float length() const {
        float len = 0.0f;
        for (size_t i = 1; i < pts.size(); ++i) len += std::sqrt(sqDist(pts[i - 1], pts[i]));
        return len;
    }
};

// Seconds to fly `cells` cells at DRONE_SPEED_MPS
inline float flightSeconds(float cells) { return cells * SIM_CELL_SIZE_M / DRONE_SPEED_MPS; }

// Length of a closed route (including the wrap from the last waypoint back to the first), cells
inline float closedRouteLength(const std::vector<Vec2>& route) {
    float len = 0.0f;
    for (size_t i = 0, n = route.size(); n > 1 && i < n; ++i) len += std::sqrt(sqDist(route[i], route[(i + 1) % n]));
    return len;
}

struct MissionPlan {
    std::vector<int> label;                        // Voronoi region per cell of the outer zone (the ACO region);
                                                   // -1 everywhere else, including the inner band
    std::vector<std::uint8_t> zone;                // Per cell: 1 inner band (sweep), 2 outer zone (ACO), 0 neither
    GridRect bounds;                               // Every cell with zone != 0 (so every labelled cell) lies inside
    std::vector<Region> regions;

    // --- Sweep ---
    std::vector<std::vector<Vec2>> perimeterLoops; // Perimeter line, as traced (tracePerimeterLine)
    std::vector<SweepLane> lanes;                  // Every lane, in sequence order along the perimeter line
    std::vector<int> owner;                        // Per cell: sector whose lane is nearest (-1 outside the band)
    std::vector<std::vector<Vec2>> sweepRoutes;    // One closed route per sector
    std::vector<int> sectorDrone;                  // Id of the drone flying each sector (-1: a new drone is needed)
    std::vector<float> sweepLapS;                  // Lap time of each sector's closed route, seconds
    int sweepNeeded = 0;                           // Sweep drones needed to keep every lap within MAX_SWEEP_LAP_S
                                                   // (clamped to MIN_SWEEP_DRONES .. MAX_SWEEP_DRONES)
    bool sweepOverLimit = false;                   // Some lap is over MAX_SWEEP_LAP_S (only at MAX_SWEEP_DRONES)
    bool sweepShrinkOk = false;                    // One drone fewer (split built and timed) would keep every
                                                   // lap within SHRINK_MARGIN of the cap on this re-plan
    bool sweepShrank = false;                      // This plan has one sector fewer than the fleet had drones:
                                                   // the drone left without a sector is retired

    // --- ACO ---
    std::vector<Vec2> acoOrder;                    // Points of interest in the order the ant colony chose (before
                                                   // the fire-safe legs are added); reused while the perimeter holds
    std::vector<Vec2> acoTour;                     // Closed tour through the outer points of interest
    bool acoKept = false;                          // Regions and tour order carried over from the previous plan
                                                   // because the perimeter had not moved (SweepFleet::perimeterUnchanged)
    float acoLoopS = 0.0f;                         // Time for ONE drone to fly the whole closed tour, seconds
    int acoNeeded = MIN_ACO_DRONES;                // ceil(acoLoopS / ACO_REVISIT_S), clamped to MIN/MAX_ACO_DRONES
    bool acoOverLimit = false;                     // Even MAX_ACO_DRONES can't keep the revisit within ACO_REVISIT_S
};
