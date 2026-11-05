#pragma once
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

    bool cellOk(int x, int y) const;
    bool lineClear(Vec2 a, Vec2 b) const;

    // Waypoints to fly from a to b (excludes a, ends with b)
    std::vector<Vec2> leg(Vec2 a, Vec2 b);

private:
    int snap(Vec2 p) const;                 // Nearest flyable cell to p
    float heuristic(int a, int b) const;    // Octile distance

    const DistanceTransform& edt_;
    std::vector<float> g_;
    std::vector<int> came_;
    std::vector<char> closed_;
};
