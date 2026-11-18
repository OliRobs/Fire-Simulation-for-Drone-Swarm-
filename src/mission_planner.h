#pragma once
#include <random>
#include <vector>
#include "distance_transform.h"
#include "mission_types.h"
#include "safe_planner.h"

// The sweep drones the plan is made for. The planner sizes the fleet (never fewer sectors than drones here:
// growing is immediate) and matches its sectors to these drones, keeping each on the land it had before.
struct SweepFleet {
    std::vector<int> ids;          // Active sweep drones (stable ids), in launch order
    std::vector<Vec2> positions;   // Where each of them is now (same order)
    const MissionPlan* previous = nullptr;   // Last plan (its owner map and sectorDrone), for stable ownership
    int shrinkStreak = 0;                    // Re-plans in a row (before this one) on which one drone fewer
                                             // would have kept every lap within SHRINK_MARGIN of the cap
};

// Builds the whole survey plan from the live fire (standoff) and the committed perimeter (zones).
// Sweep: as many sectors as needed to keep every lap within MAX_SWEEP_LAP_S, but at least as many as the fleet
// has drones, clamped to MIN_SWEEP_DRONES .. MAX_SWEEP_DRONES. plan.sectorDrone gives each sector's drone id,
// or -1 where a new drone has to be launched. The only time it plans for fewer sectors than drones is the shrink
// rule: when one fewer would keep every lap within SHRINK_MARGIN of the cap on SHRINK_REPLANS re-plans in a row
// (this one included), it uses that split, and the drone left without a sector is to be retired.
MissionPlan planMission(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt,
                        SafePlanner& planner, std::mt19937& rng, const SweepFleet& fleet);

// Same, for the starting fleet (SWEEP_DRONES sweep drones with ids 0, 1, ...), e.g. for tests
MissionPlan planMission(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt,
                        SafePlanner& planner, std::mt19937& rng);
