#pragma once
#include <random>
#include <vector>
#include "distance_transform.h"
#include "mission_types.h"
#include "safe_planner.h"

// The sweep drones the plan is made for. The planner splits the lanes into one sector per drone.
struct SweepFleet {
    std::vector<int> ids;          // Active sweep drones (stable ids), in launch order
    std::vector<Vec2> positions;   // Where each of them is now (same order)
};

// Builds the whole survey plan from the live fire (standoff) and the committed perimeter (zones).
MissionPlan planMission(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt,
                        SafePlanner& planner, std::mt19937& rng, const SweepFleet& fleet);

// Same, for the starting fleet (SWEEP_DRONES sweep drones with ids 0, 1, ...), e.g. for tests
MissionPlan planMission(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt,
                        SafePlanner& planner, std::mt19937& rng);
