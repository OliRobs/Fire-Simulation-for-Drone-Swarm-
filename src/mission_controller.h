#pragma once
#include <cstdint>
#include <random>
#include <utility>
#include <vector>
#include "distance_transform.h"
#include "drone.h"
#include "mission_types.h"
#include "safe_planner.h"

// Owns the drones and their current mission plan. Re-plans (at most every REPLAN_INTERVAL frames) when the
// fire changes, and between re-plans keeps every drone on fire-safe legs: if the fire grows across the leg a
// drone is flying, that drone is re-routed straight away (or holds position if nothing is reachable).
// The fleet follows the fire: each re-plan launches new sweep / ACO drones at once when more are needed, and
// retires one (flown home to the launch point, then removed) only after SHRINK_REPLANS re-plans in a row on
// which one fewer would have done (config.h, "Shrinking the fleet").
// Speed (fraction of DRONE_SPEED) for an ACO drone that is `gap` cells along the tour behind the drone ahead,
// when an even share of the tour is `even` cells: full speed within ACO_SPACING_TOLERANCE of an even share or
// further, otherwise in proportion to the gap, but never below ACO_SPACING_MIN_SPEED.
float acoSpacingSpeed(float gap, float even);

class MissionController {
public:
    MissionController(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt, Vec2 launchPoint);

    void markFireChanged() { planDirty_ = true; }
    void update();              // One frame: re-plan if due, then move every drone

    bool takePlanChanged();     // True once after each re-plan (renderer rebuilds its overlay)
    const MissionPlan& plan() const { return plan_; }
    const std::vector<Drone>& drones() const { return drones_; }
    int activeCount(bool aco) const;   // Active (not returning) drones of one role

    // Gaps along the ACO tour (cells) from each ACO drone that is flying the tour to the next one ahead of it,
    // in order round the tour. Drones still in transit to the tour, holding or returning are not counted.
    std::vector<float> acoGaps() const;

private:
    void replan();

    int launchDrone(bool aco);                          // New drone at the launch point; returns its id
    int chooseDroneToRetire(bool aco, int spareId) const; // Which drone leaves when the fleet shrinks
    void retire(int id);                                // Send a drone home; removed when it gets there
    Drone* findDrone(int id);                           // nullptr if there is no such drone

    // ACO drones flying the tour now, as (index into drones_, arc position along the tour in cells), sorted by
    // arc position, i.e. in the order they fly round the tour
    std::vector<std::pair<size_t, float>> acoTourPositions() const;
    void updateAcoSpacing();   // Sets each ACO drone's speedScale for this frame from its gap to the drone ahead
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
    int sweepShrinkStreak_ = 0;  // Re-plans in a row on which one sweep drone fewer would have done
    int acoShrinkStreak_ = 0;    // Re-plans in a row on which one ACO drone fewer would have done
    std::vector<float> acoArc_;  // Arc length from the ACO tour's first waypoint to each waypoint (cells)
    float acoTotal_ = 0.0f;      // Length of the closed ACO tour (cells)
    bool planned_ = false;                      // A plan has been made (plan_ is a real previous plan)
    std::uint64_t plannedPerimVersion_ = 0;     // perim_.version() when plan_ was made (see replan())
    bool planDirty_ = true;
    bool planChanged_ = false;
    int frame_ = 0;
    int lastPlanFrame_ = -REPLAN_INTERVAL;
};
