#pragma once
#include <cstdint>
#include <vector>
#include "config.h"
#include "distance_transform.h"

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
    std::vector<std::uint8_t> committed_, trigger_;
    DistanceTransform perimEdt_, triggerEdt_;
    GridRect fireSeen_;   // Union of every fire rectangle passed in so far: the committed fire and any
                          // trigger_ cells left from an earlier update all lie inside it
};
