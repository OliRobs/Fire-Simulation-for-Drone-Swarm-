#pragma once
#include <cstdint>
#include <vector>
#include "config.h"
#include "distance_transform.h"
#include "types.h"

// The perimeter (PERIMETER_DISTANCE cells out), measured from a *committed* copy of the fire. The committed fire only absorbs live
// fire near sections where the front has come within PERIMETER_UPDATE_DISTANCE of the ring, so the rest of
// the perimeter (and the survey zones measured from it) stays put.
class Perimeter {
public:
    explicit Perimeter(const std::vector<std::uint8_t>& initialFire);

    // Call when the fire has grown. Returns true if any section of the perimeter moved.
    // fireBounds must contain every non-zero cell of liveFire (FireSimulation::fireBounds()): the update
    // then only looks inside it rather than at all GRID_WIDTH x GRID_HEIGHT cells.
    bool update(const std::vector<std::uint8_t>& liveFire, const GridRect& fireBounds);
    // Same, finding the fire's bounding box by scanning liveFire (for callers without one, e.g. tests)
    bool update(const std::vector<std::uint8_t>& liveFire);

    const DistanceTransform& distance() const { return perimEdt_; } // Distance from the committed fire
    std::vector<int> ringCells() const;                              // Cells ~PERIMETER_DISTANCE out

private:
    std::vector<std::uint8_t> committed_, trigger_;   // trigger_ is all 0 between calls (see update())
    std::vector<std::uint8_t> lastLive_;   // Copy of the live-fire mask from the last call (or the starting fire),
                                           // kept up to date inside fireSeen_ (outside it, both masks are all 0)
    std::vector<int> newTriggers_;         // Scratch: cells that trigger an update in this call
    DistanceTransform perimEdt_, triggerEdt_;
    GridRect fireSeen_;   // Union of every fire rectangle passed in so far: the committed fire and any
                          // trigger_ cells left from an earlier update all lie inside it
};

// ----------------------------------------------------------------------------
// Perimeter line
// ----------------------------------------------------------------------------
// The line `level` cells out from the committed fire (by default the PERIMETER_DISTANCE ring), traced from
// the distance field as ordered closed polylines. ringCells() gives the same ring as an unordered set of
// cells; the sweep planner needs it as a line it can walk along (to seed lanes at even spacing).
//
// Each separate fire, or a large unburnt island, can give its own loop. Every loop runs the same way round
// the fire: anticlockwise as seen on screen (north up), so the fire is always on the left of the direction of
// travel. In grid coordinates (y pointing down) such a loop has a NEGATIVE signed area (loopSignedArea).
// The last point is not repeated: each loop wraps from its last point back to its first.
// Each loop starts at its topmost point (smallest y, then smallest x), so the start stays put as long as the
// fire does not change near there.
std::vector<std::vector<Vec2>> tracePerimeterLine(const DistanceTransform& perimEdt,
                                                  float level = static_cast<float>(PERIMETER_DISTANCE));

// Signed area of a closed polyline in grid coordinates (shoelace formula). Negative for a loop that runs
// anticlockwise on screen (see above).
double loopSignedArea(const std::vector<Vec2>& loop);

// Length of a closed polyline, including the closing segment back to the start
float loopLength(const std::vector<Vec2>& loop);
