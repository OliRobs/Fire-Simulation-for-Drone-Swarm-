#pragma once
#include <vector>
#include "distance_transform.h"
#include "safe_planner.h"
#include "types.h"

struct Drone {
    Vec2 pos{ 0.0f, 0.0f };
    std::vector<Vec2> route;   // Closed loop of waypoints (grid units, cell centres)
    std::vector<Vec2> transit; // Fire-safe path flown first, to reach the loop after a re-plan or re-route
    size_t target = 0;         // Index into route of the waypoint being flown to
    bool isACO = false;
    bool holding = false;      // No safe way to any nearby waypoint: hover in place until a route opens up
    int nextRetryFrame = 0;    // While holding, the frame at which re-routing is next attempted
};

size_t nearestWaypoint(const std::vector<Vec2>& route, Vec2 p);

// Fly at DRONE_SPEED: first the fire-safe transit path (if any), then round the looping route.
// If a planner is given, every segment is checked against the live fire before the drone starts along it;
// the drone stops at the waypoint instead of flying into a segment the fire has since blocked.
// Returns false if it stopped because of a blocked segment (the caller should re-route the drone).
bool stepDrone(Drone& d, const SafePlanner* planner = nullptr);

// Hard safety net: a drone closer than HARD_MIN_FIRE_DISTANCE to the fire is pushed straight back out.
void keepClearOfFire(Drone& d, const DistanceTransform& liveEdt);
