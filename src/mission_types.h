#pragma once
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>
#include "types.h"

// One Voronoi region of the search area. Regions now only supply the ACO points of interest (one in the outer
// zone of each region); the sweep is planned from lanes along the perimeter line instead (SweepLane).
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

struct MissionPlan {
    std::vector<int> label;                        // Voronoi region per cell (-1 = outside search area)
    std::vector<std::uint8_t> zone;                // Per cell: 1 inner band (sweep), 2 outer zone (ACO), 0 neither
    GridRect bounds;                               // Every cell with label >= 0 lies inside this rectangle
    std::vector<Region> regions;

    // --- Sweep ---
    std::vector<std::vector<Vec2>> perimeterLoops; // Perimeter line, as traced (tracePerimeterLine)
    std::vector<SweepLane> lanes;                  // Every lane, in sequence order along the perimeter line
    std::vector<int> owner;                        // Per cell: sector whose lane is nearest (-1 outside the band)
    std::vector<std::vector<Vec2>> sweepRoutes;    // One closed route per sector
    std::vector<int> sectorDrone;                  // Id of the drone flying each sector (-1: not yet assigned)

    // --- ACO ---
    std::vector<Vec2> acoTour;                     // Closed tour through the outer points of interest
};
