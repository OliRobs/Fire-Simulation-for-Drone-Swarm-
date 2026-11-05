#pragma once
#include <vector>
#include "distance_transform.h"
#include "types.h"

struct Drone {
    Vec2 pos{ 0.0f, 0.0f };
    std::vector<Vec2> route;   // Closed loop of waypoints (grid units, cell centres)
    std::vector<Vec2> transit; // Fire-safe path flown first, to reach the loop after a re-plan
    size_t target = 0;
    bool isACO = false;
};

size_t nearestWaypoint(const std::vector<Vec2>& route, Vec2 p);

// Fly at DRONE_SPEED: first the fire-safe transit path (if any), then round the looping route.
void stepDrone(Drone& d);

// Hard safety net: a drone closer than HARD_MIN_FIRE_DISTANCE to the fire is pushed straight back out.
void keepClearOfFire(Drone& d, const DistanceTransform& liveEdt);
