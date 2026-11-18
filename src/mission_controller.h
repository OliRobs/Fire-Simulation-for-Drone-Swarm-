#pragma once
#include <random>
#include <vector>
#include "distance_transform.h"
#include "drone.h"
#include "mission_types.h"
#include "safe_planner.h"

// Owns the drones and their current mission plan. Re-plans (at most every REPLAN_INTERVAL frames) when the
// fire changes, and between re-plans keeps every drone on fire-safe legs: if the fire grows across the leg a
// drone is flying, that drone is re-routed straight away (or holds position if nothing is reachable).
class MissionController {
public:
    MissionController(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt, Vec2 launchPoint);

    void markFireChanged() { planDirty_ = true; }
    void update();              // One frame: re-plan if due, then move every drone

    bool takePlanChanged();     // True once after each re-plan (renderer rebuilds its overlay)
    const MissionPlan& plan() const { return plan_; }
    const std::vector<Drone>& drones() const { return drones_; }
    int activeCount(bool aco) const;   // Active (not returning) drones of one role

private:
    void replan();

    int launchDrone(bool aco);                          // New drone at the launch point; returns its id
    std::vector<size_t> activeIndices(bool aco) const;  // Positions in drones_ of the active drones of a role

    // Give d a fire-safe transit path to its current route target. If that target can't be reached, try
    // the following waypoints along the route. If none can be reached the drone is put into holding.
    void routeTo(Drone& d);

    const DistanceTransform& live_;
    const DistanceTransform& perim_;
    SafePlanner planner_;
    std::mt19937 rng_{ 2024 };
    MissionPlan plan_;
    Vec2 launch_;                // Where drones launch from (and, later, return to)
    std::vector<Drone> drones_;
    int nextId_ = 0;             // Id for the next drone launched
    bool planDirty_ = true;
    bool planChanged_ = false;
    int frame_ = 0;
    int lastPlanFrame_ = -REPLAN_INTERVAL;
};
