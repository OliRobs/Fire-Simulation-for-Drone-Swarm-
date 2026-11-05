#include "mission_controller.h"
#include <algorithm>
#include <limits>
#include "mission_planner.h"

// ============================================================================
// Autonomous Multi-Drone Mission Controller
//
// Manages the operational lifecycle of a drone fleet deployed for wildfire 
// monitoring and coverage:
// 1. Manages high-level mission re-planning (combining Sector Sweeps & ACO tours).
// 2. Synchronizes path planning with real-time fire perimeter developments.
// 3. Handles path execution, transit routing, and emergency fire-avoidance safety steps.
// ============================================================================

// Constructor: Initializes grid EDT handles, planner state, and launches drones at base station
MissionController::MissionController(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt,
    Vec2 launchPoint)
    : live_(liveEdt), perim_(perimEdt), planner_(liveEdt), drones_(SWEEP_DRONES + ACO_DRONES) {
    for (int i = 0; i < SWEEP_DRONES + ACO_DRONES; ++i) {
        drones_[i].pos = launchPoint;          // Spawn all drones at base station coordinates
        drones_[i].isACO = i >= SWEEP_DRONES; // First 'SWEEP_DRONES' are sector sweepers, remainder are ACO perimeter monitors
    }
}

// ----------------------------------------------------------------------------
// Plan State Query
// ----------------------------------------------------------------------------
// Checks if the active mission plan was updated and resets the notification flag
bool MissionController::takePlanChanged() {
    const bool changed = planChanged_;
    planChanged_ = false;
    return changed;
}

// ----------------------------------------------------------------------------
// Main Controller Update (Per Frame / Tick)
// ----------------------------------------------------------------------------
// Advances mission clock, evaluates periodic re-planning needs, and updates 
// spatial positions / safety nets for all active drones.
void MissionController::update() {
    ++frame_;

    // Trigger periodic mission re-plan if fire front changed and re-plan timer expired
    if (planDirty_ && frame_ - lastPlanFrame_ >= REPLAN_INTERVAL) {
        replan();
    }

    // Step forward physical drone positions and enforce dynamic fire-avoidance safety net
    for (Drone& d : drones_) {
        stepDrone(d);
        keepClearOfFire(d, live_);
    }
}

// ----------------------------------------------------------------------------
// Mission Re-Planning Engine
// ----------------------------------------------------------------------------
// Computes updated sector coverage paths and Ant Colony Optimization (ACO) outer-tour 
// trajectories based on latest live distance fields and fire perimeter boundaries.
void MissionController::replan() {
    // Seed random generator with fixed seed to guarantee deterministic path planning between steps
    rng_.seed(2024);

    // Generate new mission plan containing sector sweep paths and an ACO perimeter loop
    plan_ = planMission(live_, perim_, planner_, rng_);
    lastPlanFrame_ = frame_;
    planDirty_ = false;
    planChanged_ = true;

    // --- 1. Assign Sector Sweep Routes ---
    // Sector-sweep drones each take responsibility for one inner-zone Voronoi sector
    for (int k = 0; k < SWEEP_DRONES; ++k) {
        drones_[k].route = plan_.sweepRoutes[k];
        drones_[k].target = nearestWaypoint(drones_[k].route, drones_[k].pos);
    }

    // --- 2. Assign ACO Outer-Perimeter Loop Routes ---
    // Distribute ACO perimeter-monitoring drones along the global optimal TSP perimeter tour
    const int n = static_cast<int>(plan_.acoTour.size());
    std::vector<char> claimed(n, 0); // Track claimed starting waypoints to prevent initial drone collisions

    for (int k = 0; k < ACO_DRONES; ++k) {
        Drone& d = drones_[SWEEP_DRONES + k];
        d.route = plan_.acoTour;
        if (n == 0) continue; // Skip if no perimeter route points available

        const bool respectClaims = n >= ACO_DRONES;
        int closest = 0;
        float bd = std::numeric_limits<float>::max();

        // Find nearest unclaimed waypoint on the ACO loop route
        for (int i = 0; i < n; ++i) {
            if (respectClaims && claimed[i]) continue;
            const float dd = sqDist(d.pos, plan_.acoTour[i]);
            if (dd < bd) {
                bd = dd;
                closest = i;
            }
        }
        claimed[closest] = 1; // Mark waypoint as claimed by drone 'k'

        // Alternate traversal direction (clockwise vs. counter-clockwise) for even fleet dispersion
        if (k % 2 == 1) {
            std::reverse(d.route.begin(), d.route.end());
            closest = n - 1 - closest; // Invert entry target index for reversed route vector
        }
        d.target = static_cast<size_t>(closest);
    }

    // --- 3. Generate Fire-Safe Transit Paths ---
    // Calculate safe A* transit paths connecting each drone's current position to its assigned starting waypoint
    for (Drone& d : drones_) {
        d.transit.clear();
        if (!d.route.empty()) {
            d.transit = planner_.leg(d.pos, d.route[d.target]);
        }
    }
}