#pragma once
#include <random>
#include "distance_transform.h"
#include "mission_types.h"
#include "safe_planner.h"

// Builds the whole survey plan from the live fire (standoff) and the committed perimeter (zones).
MissionPlan planMission(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt,
                        SafePlanner& planner, std::mt19937& rng);
