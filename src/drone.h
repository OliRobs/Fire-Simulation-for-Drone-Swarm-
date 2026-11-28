#pragma once
#include <vector>
#include "distance_transform.h"
#include "safe_planner.h"
#include "types.h"

// What a drone is doing. Active drones fly their assigned route; a Returning drone has been retired (the fleet
// shrank) and is flying a fire-safe path back to the launch point, where the controller removes it.
enum class DroneState { Active, Returning };

struct Drone {
    int id = -1;               // Stable identity, unique for the whole run (the fleet grows and shrinks, so a
                               // drone's position in the controller's list is not a stable way to refer to it)
    DroneState state = DroneState::Active;
    Vec2 pos{ 0.0f, 0.0f };
    std::vector<Vec2> route;   // Closed loop of waypoints (grid units, cell centres)
    std::vector<Vec2> transit; // Fire-safe path flown first, to reach the loop after a re-plan or re-route
    size_t target = 0;         // Index into route of the waypoint being flown to
    bool isACO = false;
    bool holding = false;      // No safe way to any nearby waypoint: hover in place until a route opens up
    float speedScale = 1.0f;   // Fraction of DRONE_SPEED flown this frame. Below 1 only for an ACO drone that is
                               // too close behind the next one, to open the gap back up (MissionController)
    int nextRetryFrame = 0;    // While holding, the frame at which re-routing is next attempted
};

size_t nearestWaypoint(const std::vector<Vec2>& route, Vec2 p);

// Fly at DRONE_SPEED x d.speedScale: first the fire-safe transit path (if any), then round the looping route.
// If a planner is given, every segment is checked against the live fire before the drone starts along it;
// the drone stops at the waypoint instead of flying into a segment the fire has since blocked.
// Returns false if it stopped because of a blocked segment (the caller should re-route the drone).
bool stepDrone(Drone& d, const SafePlanner* planner = nullptr);

// Hard safety net: a drone closer than HARD_MIN_FIRE_DISTANCE to the fire is pushed straight back out
// (kept on liveEdt's grid).
void keepClearOfFire(Drone& d, const DistanceTransform& liveEdt);
