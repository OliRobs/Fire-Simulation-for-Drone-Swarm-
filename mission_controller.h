#pragma once
#include <random>
#include <vector>
#include "distance_transform.h"
#include "drone.h"
#include "mission_types.h"
#include "safe_planner.h"

// Owns the drones and the current survey plan. Re-plans (rate-limited) when the fire has changed, hands the
// new routes to the drones, and flies them each frame.
class MissionController {
public:
    MissionController(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt, Vec2 launchPoint);

    void markFireChanged() { planDirty_ = true; }
    void update();              // One frame: re-plan if due, then move every drone

    bool takePlanChanged();     // True once after each re-plan (renderer rebuilds its overlay)
    const MissionPlan& plan() const { return plan_; }
    const std::vector<Drone>& drones() const { return drones_; }

private:
    void replan();

    const DistanceTransform& live_;
    const DistanceTransform& perim_;
    SafePlanner planner_;
    std::mt19937 rng_{ 2024 };
    MissionPlan plan_;
    std::vector<Drone> drones_;
    bool planDirty_ = true;
    bool planChanged_ = false;
    int frame_ = 0;
    int lastPlanFrame_ = -REPLAN_INTERVAL;
};
