#pragma once
#include <cstdint>
#include <vector>
#include "distance_transform.h"
#include "types.h"

// Fire-safe path finding. A cell is flyable if it is at least HARD_MIN_FIRE_DISTANCE from the fire, so legs
// may cut through the wider MIN_FIRE_DISTANCE survey standoff whenever the path stays clear of the fire itself.
// leg() returns the straight line when it is clear, otherwise an A* path around the fire pulled taut with
// line-of-sight checks.
class SafePlanner {
public:
    explicit SafePlanner(const DistanceTransform& liveEdt);

    // True if cell (x, y) is on the map and outside the hard no-fly zone
    bool cellOk(int x, int y) const;

    // True if the straight segment a -> b never enters the hard no-fly zone
    bool lineClear(Vec2 a, Vec2 b) const;

    // Fire-safe waypoints to fly from a to b. The start point a is not included.
    //  - Normally the path ends exactly at b.
    //  - If b itself is inside the hard no-fly zone, the path ends at the nearest flyable cell to b.
    //  - If there is no safe path at all (b is cut off by fire), the result is EMPTY.
    //    Callers must treat empty as "do not fly there"; it is never replaced by a straight line.
    std::vector<Vec2> leg(Vec2 a, Vec2 b);

private:
    int snap(Vec2 p) const;                 // Nearest flyable cell to p, or -1 if none within 40 cells
    float heuristic(int a, int b) const;    // Octile distance (admissible for 8-connected moves)

    const DistanceTransform& edt_;

    // A* working storage, reused between calls. Instead of clearing all 250k entries on every search,
    // each search bumps search_ and an entry only counts as valid if its stamp equals the current search.
    std::vector<float> g_;                  // Cost from the start cell
    std::vector<int> came_;                 // Parent cell on the best path found so far
    std::vector<std::uint32_t> openStamp_;  // == search_ when g_/came_ hold values for this search
    std::vector<std::uint32_t> closedStamp_;// == search_ once the cell has been expanded in this search
    std::uint32_t search_ = 0;
};
