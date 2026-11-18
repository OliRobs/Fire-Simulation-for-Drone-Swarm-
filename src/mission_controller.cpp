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

// Constructor: Initializes grid EDT handles, planner state, and launches the starting fleet at the base
// station: SWEEP_DRONES sweep drones and ACO_DRONES ACO drones, each with its own stable id.
MissionController::MissionController(const DistanceTransform& liveEdt, const DistanceTransform& perimEdt,
    Vec2 launchPoint)
    : live_(liveEdt), perim_(perimEdt), planner_(liveEdt), launch_(launchPoint) {
    for (int i = 0; i < SWEEP_DRONES; ++i) launchDrone(false);
    for (int i = 0; i < ACO_DRONES; ++i) launchDrone(true);
}

// ----------------------------------------------------------------------------
// Launch One Drone
// ----------------------------------------------------------------------------
// Adds a new active drone of the given role at the launch point, with the next unused id. It has no route
// yet: the caller gives it one. Returns its id (not a reference, which the next launch could invalidate).
int MissionController::launchDrone(bool aco) {
    Drone d;
    d.id = nextId_++;
    d.isACO = aco;
    d.pos = launch_;
    drones_.push_back(d);
    return d.id;
}

// ----------------------------------------------------------------------------
// Drones by Role
// ----------------------------------------------------------------------------
// Indices into drones_ of the ACTIVE drones of one role (sweep or ACO), in launch order. Planning code goes
// through this instead of assuming the first SWEEP_DRONES entries are the sweepers.
std::vector<size_t> MissionController::activeIndices(bool aco) const {
    std::vector<size_t> out;
    for (size_t i = 0; i < drones_.size(); ++i) {
        if (drones_[i].isACO == aco && drones_[i].state == DroneState::Active) out.push_back(i);
    }
    return out;
}

int MissionController::activeCount(bool aco) const { return static_cast<int>(activeIndices(aco).size()); }

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

        // A returning drone whose path home ended short of the launch point (the fire was in the way when it
        // was planned) asks for a new one
        if (d.state == DroneState::Returning && !d.holding && d.transit.empty() && sqDist(d.pos, launch_) > 1e-4f) {
            routeTo(d);
        }
    }

    // Retired drones that have got back to the launch point leave the fleet
    drones_.erase(std::remove_if(drones_.begin(), drones_.end(), [&](const Drone& d) {
        return d.state == DroneState::Returning && d.transit.empty() && sqDist(d.pos, launch_) <= 1e-4f;
    }), drones_.end());
}

// ----------------------------------------------------------------------------
// Re-route One Drone
// ----------------------------------------------------------------------------
// Finds a fire-safe transit path from the drone's position to its route target. Waypoints that can only be
// reached by crossing the fire are skipped (trying up to MAX_REROUTE_TRIES of them in route order).
// For a returning (retired) drone the target is the launch point.
void MissionController::routeTo(Drone& d) {
    d.transit.clear();

    // A retired drone flies home: a fire-safe path to the launch point, or hold (and retry) if there is none
    if (d.state == DroneState::Returning) {
        std::vector<Vec2> path = planner_.leg(d.pos, launch_);
        d.holding = path.empty() && sqDist(d.pos, launch_) > 1e-4f;
        d.nextRetryFrame = frame_ + HOLD_RETRY_FRAMES;
        d.transit = std::move(path);
        return;
    }

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

    // The sweep drones the lanes are split between, and the last plan (for stable ownership)
    const std::vector<size_t> sweep = activeIndices(false);
    SweepFleet fleet;
    for (size_t i : sweep) {
        fleet.ids.push_back(drones_[i].id);
        fleet.positions.push_back(drones_[i].pos);
    }
    fleet.previous = &plan_;
    fleet.shrinkStreak = sweepShrinkStreak_;

    // (The new plan is built while the old one is still in plan_, then replaces it)
    plan_ = planMission(live_, perim_, planner_, rng_, fleet);
    lastPlanFrame_ = frame_;
    planDirty_ = false;
    planChanged_ = true;

    // --- Sweep shrink counter: re-plans in a row on which one drone fewer would have done (the planner applies
    // the shrink itself on the SHRINK_REPLANS-th, and then the count starts again) ---
    sweepShrinkStreak_ = (plan_.sweepShrinkOk && !plan_.sweepShrank) ? sweepShrinkStreak_ + 1 : 0;

    // --- Sweep fleet growth: a sector no existing drone was matched to gets a new drone from the launch point.
    // Growing takes effect at once (shrinking waits; see the shrink rule). ---
    for (int& id : plan_.sectorDrone) {
        if (id < 0) id = launchDrone(false);
    }

    // --- Sweep fleet shrink: an active sweep drone left without a sector is retired ---
    std::vector<int> spareSweep;
    for (size_t i : activeIndices(false)) {
        if (std::find(plan_.sectorDrone.begin(), plan_.sectorDrone.end(), drones_[i].id) == plan_.sectorDrone.end())
            spareSweep.push_back(drones_[i].id);
    }
    for (int spare : spareSweep) {
        const int out = chooseDroneToRetire(false, spare);
        // If the choice is a drone that does have a sector, the spare drone takes that sector over
        for (int& id : plan_.sectorDrone) if (id == out) id = spare;
        retire(out);
    }

    // --- Sweep drones: each flies the sector the plan gave its id, joined at the waypoint nearest it. A drone
    // left without a sector (fewer lanes than drones) waits in place with no route. ---
    for (size_t i : activeIndices(false)) {
        Drone& d = drones_[i];
        const auto it = std::find(plan_.sectorDrone.begin(), plan_.sectorDrone.end(), d.id);
        if (it == plan_.sectorDrone.end()) { d.route.clear(); continue; }
        d.route = plan_.sweepRoutes[it - plan_.sectorDrone.begin()];
        d.target = nearestWaypoint(d.route, d.pos);
    }

    // --- ACO drones: all share the same closed tour, flown in the same direction, evenly spaced ---
    // Spacing them by tour length (rather than each grabbing its nearest waypoint) stops them flying in a
    // bunch, so any point on the tour is revisited every (tour length / ACO count) instead of in bursts.
    // ACO fleet size. Growth: launch new drones at once from the launch point if the tour needs more for the
    // revisit limit (plan.acoNeeded). Shrink: one drone fewer once the remaining drones would have kept the
    // revisit time within SHRINK_MARGIN of the limit on SHRINK_REPLANS re-plans in a row. Slots are then worked
    // out for the new count, so the drones are re-spaced evenly.
    int numAco = activeCount(true);
    if (plan_.acoNeeded > numAco) {
        for (int k = numAco; k < plan_.acoNeeded; ++k) launchDrone(true);
        numAco = plan_.acoNeeded;
        acoShrinkStreak_ = 0;
    }
    else if (numAco > MIN_ACO_DRONES && plan_.acoLoopS / (numAco - 1) <= SHRINK_MARGIN * ACO_REVISIT_S) {
        if (++acoShrinkStreak_ >= SHRINK_REPLANS) {
            --numAco;              // The drone left without a slot below is retired
            acoShrinkStreak_ = 0;
        }
    }
    else {
        acoShrinkStreak_ = 0;
    }

    const std::vector<size_t> aco = activeIndices(true);
    const std::vector<Vec2>& tour = plan_.acoTour;
    const int n = static_cast<int>(tour.size());
    std::vector<int> spareAco;   // ACO drones left without a slot (fleet shrinking)
    if (n > 0 && numAco > 0) {
        // Arc length from waypoint 0 to each waypoint, and the total closed-loop length
        std::vector<float> arc(n, 0.0f);
        for (int i = 1; i < n; ++i) arc[i] = arc[i - 1] + std::sqrt(sqDist(tour[i - 1], tour[i]));
        const float total = arc[n - 1] + std::sqrt(sqDist(tour[n - 1], tour[0]));

        // Anchor the spacing at the waypoint closest to the ACO drones' average position, so the group
        // as a whole joins the tour at its closest point
        Vec2 mean(0.0f, 0.0f);
        for (size_t i : aco) mean += drones_[i].pos;
        mean = mean / static_cast<float>(aco.size());
        const int anchor = static_cast<int>(nearestWaypoint(tour, mean));

        // Slot j is the waypoint closest to (anchor + j/numAco of the way round the tour)
        std::vector<int> slots(numAco);
        for (int j = 0; j < numAco; ++j) {
            const float want = std::fmod(arc[anchor] + total * j / numAco, std::max(total, 1e-6f));
            int best = 0;
            float bd = std::numeric_limits<float>::max();
            for (int i = 0; i < n; ++i) {
                const float dd = std::fabs(arc[i] - want);
                if (dd < bd) { bd = dd; best = i; }
            }
            slots[j] = best;
        }

        // Each drone takes the closest slot not yet claimed by another drone. When the fleet is shrinking there
        // is one slot fewer than drones, and the drone left without one is retired.
        std::vector<char> claimed(numAco, 0);
        for (size_t i : aco) {
            Drone& d = drones_[i];
            int pick = -1;
            float bd = std::numeric_limits<float>::max();
            for (int j = 0; j < numAco; ++j) {
                if (claimed[j]) continue;
                const float dd = sqDist(d.pos, tour[slots[j]]);
                if (dd < bd) { bd = dd; pick = j; }
            }
            if (pick < 0) { spareAco.push_back(d.id); continue; }
            claimed[pick] = 1;
            d.route = tour;
            d.target = static_cast<size_t>(slots[pick]);
        }
    }
    else {
        // No tour: no slots, but the count still follows the rule (the last drones in launch order are spare)
        for (size_t k = 0; k < aco.size(); ++k) {
            drones_[aco[k]].route.clear();
            if (static_cast<int>(k) >= numAco) spareAco.push_back(drones_[aco[k]].id);
        }
    }
    for (int spare : spareAco) {
        const int out = chooseDroneToRetire(true, spare);
        if (out != spare) {
            // The chosen drone hands its slot over to the spare one before it leaves
            Drone* from = findDrone(out);
            Drone* to = findDrone(spare);
            if (from && to) { to->route = from->route; to->target = from->target; }
        }
        retire(out);
    }

    // --- Every drone: fire-safe transit from where it is now onto its route (or, for a retired drone, home
    // to the launch point, re-planned round the new fire) ---
    for (Drone& d : drones_) {
        d.holding = false;
        d.transit.clear();
        if (!d.route.empty() || d.state == DroneState::Returning) routeTo(d);
    }
}

// ----------------------------------------------------------------------------
// Which Drone to Retire
// ----------------------------------------------------------------------------
// Called once for each drone the fleet of one role has to lose, with the drone the re-plan left without a
// sector (sweep) or slot (ACO). Kept as its own function so the choice can later take battery level into
// account (e.g. retire the drone with the least charge left). It may return any active drone of that role:
// the caller then gives the spare drone the chosen drone's job. For now it retires the spare drone itself.
int MissionController::chooseDroneToRetire(bool /*aco*/, int spareId) const {
    return spareId;
}

// ----------------------------------------------------------------------------
// Retire One Drone
// ----------------------------------------------------------------------------
// The drone stops surveying and flies a fire-safe route back to the launch point, where update() removes it
// from the fleet. If there is no safe route it holds in place (keepClearOfFire still applies) and retries.
void MissionController::retire(int id) {
    Drone* d = findDrone(id);
    if (!d) return;
    d->state = DroneState::Returning;
    d->route.clear();
    d->target = 0;
    d->holding = false;
    routeTo(*d);
}

Drone* MissionController::findDrone(int id) {
    for (Drone& d : drones_) if (d.id == id) return &d;
    return nullptr;
}
