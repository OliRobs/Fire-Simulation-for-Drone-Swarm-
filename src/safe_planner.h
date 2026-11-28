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
    // Plans against liveEdt, on liveEdt's grid (DistanceTransform::grid())
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
    void refreshFlyable();                  // Brings fly_ up to date with the live distance field (see below)

    const DistanceTransform& edt_;
    const int W_, H_;                       // Grid size, cells (from edt_.grid(); kept here for the hot loops)

    // A* working storage, reused between calls. Instead of clearing all 250k entries on every search,
    // each search bumps search_ and an entry only counts as valid if its stamp equals the current search.
    // All of a cell's A* state in one 16-byte record, so the search touches one cache line per cell rather than
    // four separate whole-grid arrays (it jumps around the grid, so each array access was usually a cache miss).
    struct Node {
        float g;                    // Cost from the start cell
        int came;                   // Parent cell on the best path found so far
        std::uint32_t openStamp;    // == search_ when g/came hold values for this search
        std::uint32_t closedStamp;  // == search_ once the cell has been expanded in this search
    };
    std::vector<Node> nodes_;
    std::uint32_t search_ = 0;

    // A* open list, bucketed by f-cost: bucket k holds the entries with k <= f < k + 1, each bucket a small binary
    // heap of (f, cell). The smallest (f, cell) overall is the smallest in the lowest non-empty bucket, so entries
    // come out in exactly the order of a single priority queue of (f, cell) pairs (same tie-break on the cell
    // index), and A* finds exactly the same paths; but each push and pop works on a heap of a few dozen entries
    // instead of thousands. The buckets are kept between searches to reuse their memory.
    std::vector<std::vector<std::pair<float, int>>> buckets_;

    // Flyable-cell bitmap for A*: one byte per cell, 1 where cellOk() is true, on a grid padded with a border of
    // 0 (not flyable) cells all round, so neighbour tests need no bounds checks. Cell (x, y) is at
    // (y + 1) * flyW_ + (x + 1). It is only brought up to date when an A* search is about to run and the live
    // distance field has been recomputed since (DistanceTransform::version()), and only the field's window is
    // rewritten (outside it every cell is far from the fire, i.e. flyable).
    const int flyW_;                        // Padded row length, W_ + 2
    std::vector<std::uint8_t> fly_;
    bool flyBuilt_ = false;
    std::uint64_t flyVersion_ = 0;          // edt_.version() that fly_ matches
    GridRect flyWin_;                       // Window of the field fly_ was last built from
};
