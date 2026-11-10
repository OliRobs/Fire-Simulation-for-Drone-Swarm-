#include "mission_controller.h"
#include <algorithm>
#include <cmath>
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

// Config for re-routing between re-plans
namespace {
    constexpr size_t MAX_REROUTE_TRIES = 8;  // Waypoints tried (in route order) when re-routing one drone
    constexpr int HOLD_RETRY_FRAMES = 15;    // How often a holding drone looks for a safe route again
}

// ----------------------------------------------------------------------------
// Main Controller Update (Per Frame / Tick)
// ----------------------------------------------------------------------------
// Advances the mission clock, re-plans if due, then moves every drone. Between re-plans the fire keeps
// growing, so each drone's next segment is checked against the live fire as it flies (inside stepDrone);
// a blocked drone is re-routed immediately instead of being left to bounce off the no-fly zone.
void MissionController::update() {
    ++frame_;

    // Trigger periodic mission re-plan if fire front changed and re-plan timer expired
    if (planDirty_ && frame_ - lastPlanFrame_ >= REPLAN_INTERVAL) {
        replan();
    }

    for (Drone& d : drones_) {
        // Holding drones periodically retry, in case the fire has burnt out or a new plan opened a path
        if (d.holding && frame_ >= d.nextRetryFrame) routeTo(d);

        // Fly; if the next segment has become blocked, re-route this drone and fly the new path
        if (!stepDrone(d, &planner_)) {
            planDirty_ = true;   // The plan is out of date: ask for a fresh one at the next opportunity
            routeTo(d);
            stepDrone(d, &planner_);
        }

        // Last-resort safety net (the fire can also grow onto a drone that is standing still)
        keepClearOfFire(d, live_);
    }
}

// ----------------------------------------------------------------------------
// Re-route One Drone
// ----------------------------------------------------------------------------
// Finds a fire-safe transit path from the drone's position to its route target. Waypoints that can only be
// reached by crossing the fire are skipped (trying up to MAX_REROUTE_TRIES of them in route order).
void MissionController::routeTo(Drone& d) {
    d.transit.clear();
    const size_t n = d.route.size();

    for (size_t k = 0; k < std::min(n, MAX_REROUTE_TRIES); ++k) {
        const size_t idx = (d.target + k) % n;
        std::vector<Vec2> path = planner_.leg(d.pos, d.route[idx]);

        // Accept only paths that actually arrive at the waypoint. leg() may stop short at the nearest flyable
        // cell when the waypoint itself is now inside the no-fly zone; flying on from there would be blocked.
        if (!path.empty() && sqDist(path.back(), d.route[idx]) < 1e-6f) {
            d.target = idx;
            d.transit = std::move(path);
            d.holding = false;
            return;
        }
    }

    // Nothing reachable: hover where we are (keepClearOfFire still applies) and retry shortly
    d.holding = n > 0;
    d.nextRetryFrame = frame_ + HOLD_RETRY_FRAMES;
}

// ----------------------------------------------------------------------------
// Mission Re-Planning Engine
// ----------------------------------------------------------------------------
// Computes updated sector sweep routes and the ACO outer tour from the latest live fire and perimeter,
// then hands routes to the drones.
void MissionController::replan() {
    // Re-seed so the same fire state always gives the same plan (stable, reproducible routes)
    rng_.seed(2024);

    plan_ = planMission(live_, perim_, planner_, rng_);
    lastPlanFrame_ = frame_;
    planDirty_ = false;
    planChanged_ = true;

    // --- Sweep drones: one boustrophedon route each, joined at the waypoint nearest the drone ---
    for (int k = 0; k < SWEEP_DRONES; ++k) {
        drones_[k].route = plan_.sweepRoutes[k];
        drones_[k].target = nearestWaypoint(drones_[k].route, drones_[k].pos);
    }

    // --- ACO drones: all share the same closed tour, flown in the same direction, evenly spaced ---
    // Spacing them by tour length (rather than each grabbing its nearest waypoint) stops them flying in a
    // bunch, so any point on the tour is revisited every (tour length / ACO_DRONES) instead of in bursts.
    const std::vector<Vec2>& tour = plan_.acoTour;
    const int n = static_cast<int>(tour.size());
    if (n > 0) {
        // Arc length from waypoint 0 to each waypoint, and the total closed-loop length
        std::vector<float> arc(n, 0.0f);
        for (int i = 1; i < n; ++i) arc[i] = arc[i - 1] + std::sqrt(sqDist(tour[i - 1], tour[i]));
        const float total = arc[n - 1] + std::sqrt(sqDist(tour[n - 1], tour[0]));

        // Anchor the spacing at the waypoint closest to the ACO drones' average position, so the group
        // as a whole joins the tour at its closest point
        Vec2 mean(0.0f, 0.0f);
        for (int k = 0; k < ACO_DRONES; ++k) mean += drones_[SWEEP_DRONES + k].pos;
        mean = mean / static_cast<float>(ACO_DRONES);
        const int anchor = static_cast<int>(nearestWaypoint(tour, mean));

        // Slot j is the waypoint closest to (anchor + j/ACO_DRONES of the way round the tour)
        std::vector<int> slots(ACO_DRONES);
        for (int j = 0; j < ACO_DRONES; ++j) {
            const float want = std::fmod(arc[anchor] + total * j / ACO_DRONES, std::max(total, 1e-6f));
            int best = 0;
            float bd = std::numeric_limits<float>::max();
            for (int i = 0; i < n; ++i) {
                const float dd = std::fabs(arc[i] - want);
                if (dd < bd) { bd = dd; best = i; }
            }
            slots[j] = best;
        }

        // Each drone takes the closest slot not yet claimed by another drone
        std::vector<char> claimed(ACO_DRONES, 0);
        for (int k = 0; k < ACO_DRONES; ++k) {
            Drone& d = drones_[SWEEP_DRONES + k];
            d.route = tour;
            int pick = 0;
            float bd = std::numeric_limits<float>::max();
            for (int j = 0; j < ACO_DRONES; ++j) {
                if (claimed[j]) continue;
                const float dd = sqDist(d.pos, tour[slots[j]]);
                if (dd < bd) { bd = dd; pick = j; }
            }
            claimed[pick] = 1;
            d.target = static_cast<size_t>(slots[pick]);
        }
    }
    else {
        for (int k = 0; k < ACO_DRONES; ++k) drones_[SWEEP_DRONES + k].route.clear();
    }

    // --- Every drone: fire-safe transit from where it is now onto its route ---
    for (Drone& d : drones_) {
        d.holding = false;
        d.transit.clear();
        if (!d.route.empty()) routeTo(d);
    }
}
