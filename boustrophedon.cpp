#include "boustrophedon.h"
#include <algorithm>
#include <utility>

// ============================================================================
// Boustrophedon ("As the Ox Ploughs") Coverage Path Planning (CPP)
//
// Generates a lawnmower-style back-and-forth sweep pattern over a specific 
// region (e.g., a Voronoi region assigned to an autonomous agent/drone).
// 
// Characteristics:
// - Parallel vertical lanes spaced apart by `LANE_SPACING` (sensor swath width).
// - Alternating sweep directions (upwards then downwards) for efficient path continuity.
// - Handles disjoint cell runs per lane (caused by terrain obstacles, clipping, or 
//   distance field boundaries) by splitting each lane into individual contiguous segments.
// ============================================================================
void appendBoustrophedon(int r, const Region& rg, const std::vector<int>& label,
    const std::vector<std::uint8_t>& zone, std::vector<Vec2>& route) {

    // ------------------------------------------------------------------------
    // Step 1: Calculate X-coordinates for vertical sweep lanes
    // ------------------------------------------------------------------------
    std::vector<int> lanes;

    // Space lanes evenly starting half a swath distance from the region's left boundary
    for (int x = rg.minX + LANE_SPACING / 2; x <= rg.maxX; x += LANE_SPACING) {
        lanes.push_back(x);
    }

    // Ensure the rightmost boundary is covered if the spacing step didn't reach it fully
    if (lanes.empty() || lanes.back() + LANE_SPACING / 2 < rg.maxX) {
        lanes.push_back(rg.maxX);
    }

    const size_t startSize = route.size();
    bool down = true; // Flag controlling current lane direction (true = top-to-bottom, false = bottom-to-top)

    // ------------------------------------------------------------------------
    // Step 2: Extract vertical runs and construct waypoints for each lane
    // ------------------------------------------------------------------------
    for (int x : lanes) {
        // Holds continuous [startY, endY] intervals of valid region cells along the current vertical line `x`
        std::vector<std::pair<int, int>> runs;
        int runStart = -1;

        // Scan vertically through the region bounding box (including 1 cell padding to close open runs)
        for (int y = rg.minY; y <= rg.maxY + 1; ++y) {
            // Check if cell (x, y) belongs to region `r` and is within the target inner zone (zone == 1)
            const bool in = y <= rg.maxY && label[y * GRID_WIDTH + x] == r && zone[y * GRID_WIDTH + x] == 1;

            if (in && runStart < 0) {
                runStart = y; // Start of a new contiguous run
            }
            else if (!in && runStart >= 0) {
                runs.emplace_back(runStart, y - 1); // End of current contiguous run
                runStart = -1;
            }
        }

        // If no valid cells were found on this vertical lane, skip to next lane
        if (runs.empty()) continue;

        // Alternate traversal order of runs depending on sweep direction
        if (!down) std::reverse(runs.begin(), runs.end());

        // Convert Y-intervals into explicit 2D floating-point waypoints centered in grid cells (+0.5f)
        for (const auto& run : runs) {
            float y0 = run.first + 0.5f, y1 = run.second + 0.5f;

            // Reverse waypoint order for bottom-to-top sweep passes
            if (!down) std::swap(y0, y1);

            route.emplace_back(x + 0.5f, y0); // Segment start point
            if (y1 != y0) route.emplace_back(x + 0.5f, y1); // Segment end point (if multi-cell run)
        }

        // Flip direction for the next lane to create the serpentine pattern
        down = !down;
    }

    // ------------------------------------------------------------------------
    // Step 3: Fallback for tiny/sliver regions
    // ------------------------------------------------------------------------
    // If no path points were appended (e.g., region too narrow for lane placement),
    // fallback to visiting the region's inner centroid directly.
    if (route.size() == startSize) {
        route.push_back(rg.innerCentroid);
    }
}