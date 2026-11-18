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

    // ACO drones too close behind the drone ahead fly slower this frame, so the spacing evens out
    updateAcoSpacing();

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
    // revisit time within SHRINK_MARGIN of the limit on SHRINK_REPLANS re-plans in a row. Drones already on the
    // tour keep their progress; new drones fill the largest gaps between them.
    int numAco = activeCount(true);
    if (plan_.acoNeeded > numAco) {
        for (int k = numAco; k < plan_.acoNeeded; ++k) launchDrone(true);
        numAco = plan_.acoNeeded;
        acoShrinkStreak_ = 0;
    }
    else if (numAco > MIN_ACO_DRONES && plan_.acoLoopS / (numAco - 1) <= SHRINK_MARGIN * ACO_REVISIT_S) {
        if (++acoShrinkStreak_ >= SHRINK_REPLANS) {
            --numAco;              // The drone in the most crowded spot below is retired
            acoShrinkStreak_ = 0;
        }
    }
    else {
        acoShrinkStreak_ = 0;
    }

    const std::vector<size_t> aco = activeIndices(true);
    const std::vector<Vec2>& tour = plan_.acoTour;
    const int n = static_cast<int>(tour.size());
    std::vector<int> spareAco;   // ACO drones to retire (fleet shrinking)
    if (n > 0 && numAco > 0) {
        // Arc length from waypoint 0 to each waypoint, and the total closed-loop length
        std::vector<float> arc(n, 0.0f);
        for (int i = 1; i < n; ++i) arc[i] = arc[i - 1] + std::sqrt(sqDist(tour[i - 1], tour[i]));
        const float total = std::max(arc[n - 1] + std::sqrt(sqDist(tour[n - 1], tour[0])), 1e-6f);
        acoArc_ = arc;          // Kept for the spacing control between re-plans
        acoTotal_ = total;

        // Waypoint whose arc position is closest to `want` (taken round the loop)
        auto waypointAt = [&](float want) {
            want = std::fmod(std::fmod(want, total) + total, total);
            int best = 0;
            float bd = std::numeric_limits<float>::max();
            for (int i = 0; i < n; ++i) {
                const float dd = std::min(std::fabs(arc[i] - want), total - std::fabs(arc[i] - want));
                if (dd < bd) { bd = dd; best = i; }
            }
            return best;
        };

        // --- Drones already flying the tour keep their progress ---
        // A drone that had a route before this re-plan locks back on at the tour node nearest to it and carries on
        // the same way round (the tour is always anticlockwise, see planMission). If it has already passed that
        // node (it is beyond it, along the leg to the next node), it goes on to the next node instead, so it
        // never turns back.
        std::vector<size_t> onTour, joining;
        for (size_t i : aco) (drones_[i].route.empty() ? joining : onTour).push_back(i);
        for (size_t i : onTour) {
            Drone& d = drones_[i];
            const size_t near = nearestWaypoint(tour, d.pos);
            const size_t next = (near + 1) % n;
            const Vec2 ahead = tour[next] - tour[near], fromNode = d.pos - tour[near];
            const bool passed = n > 1 && ahead.x * fromNode.x + ahead.y * fromNode.y > 0.0f;
            d.route = tour;
            d.target = passed ? next : near;
        }

        // --- Shrinking: retire the drone in the most crowded spot ---
        // Its two neighbours' gaps merge, so removing the drone with the shortest gap behind plus ahead of it
        // leaves the remaining drones as evenly spaced as possible. (chooseDroneToRetire may still pick another.)
        while (static_cast<int>(onTour.size() + joining.size()) > numAco) {
            if (!joining.empty()) { spareAco.push_back(drones_[joining.back()].id); joining.pop_back(); continue; }
            std::vector<size_t> order = onTour;
            std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                return arc[drones_[a].target] < arc[drones_[b].target];
            });
            const size_t m = order.size();
            size_t crowded = 0;
            float smallest = std::numeric_limits<float>::max();
            for (size_t k = 0; k < m; ++k) {
                const float prev = arc[drones_[order[(k + m - 1) % m]].target];
                const float next = arc[drones_[order[(k + 1) % m]].target];
                const float span = m == 1 ? total : std::fmod(next - prev + 2.0f * total, total);
                if (span < smallest) { smallest = span; crowded = k; }
            }
            spareAco.push_back(drones_[order[crowded]].id);
            onTour.erase(std::find(onTour.begin(), onTour.end(), order[crowded]));
        }

        // --- Drones joining the tour (just launched) ---
        if (onTour.empty()) {
            // Nobody on the tour yet (first plan): space them evenly, anchored at the node nearest the group's
            // average position, each taking the closest free slot
            Vec2 mean(0.0f, 0.0f);
            for (size_t i : joining) mean += drones_[i].pos;
            mean = mean / static_cast<float>(std::max<size_t>(joining.size(), 1));
            const float anchorArc = arc[nearestWaypoint(tour, mean)];
            const int count = static_cast<int>(joining.size());
            std::vector<char> claimed(count, 0);
            for (size_t i : joining) {
                Drone& d = drones_[i];
                int pick = 0;
                float bd = std::numeric_limits<float>::max();
                for (int j = 0; j < count; ++j) {
                    if (claimed[j]) continue;
                    const float dd = sqDist(d.pos, tour[waypointAt(anchorArc + total * j / count)]);
                    if (dd < bd) { bd = dd; pick = j; }
                }
                claimed[pick] = 1;
                d.route = tour;
                d.target = static_cast<size_t>(waypointAt(anchorArc + total * pick / count));
            }
        }
        else {
            // Others are already flying: each new drone takes the middle of the largest gap between drones, so
            // the drones already on the tour don't have to move to make room
            std::vector<float> taken;
            for (size_t i : onTour) taken.push_back(arc[drones_[i].target]);
            for (size_t i : joining) {
                std::sort(taken.begin(), taken.end());
                float bestGap = -1.0f, mid = 0.0f;
                for (size_t k = 0; k < taken.size(); ++k) {
                    const float from = taken[k];
                    const float gap = k + 1 < taken.size() ? taken[k + 1] - from : taken[0] + total - from;
                    if (gap > bestGap) { bestGap = gap; mid = from + gap * 0.5f; }
                }
                Drone& d = drones_[i];
                d.route = tour;
                d.target = static_cast<size_t>(waypointAt(mid));
                taken.push_back(arc[d.target]);
            }
        }
    }
    else {
        acoArc_.clear();
        acoTotal_ = 0.0f;
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

// ============================================================================
// ACO Spacing
//
// After a re-plan each ACO drone carries on from where it was (it never turns back), and new drones join in the
// largest gaps, so the drones can end up unevenly spaced round the new tour: several could then pass the same
// node within a short time and leave other parts unvisited for longer. Drones can't fly faster than DRONE_SPEED,
// so the spacing is restored by slowing down: a drone closer than an even share of the tour behind the drone
// ahead flies slower until the gap opens up. As the gaps add up to the whole tour, the drones with too-small
// gaps opening up is the same as the too-large gaps closing, and the group settles evenly spaced, all at full
// speed again.
// ============================================================================

float acoSpacingSpeed(float gap, float even) {
    if (even <= 0.0f || gap >= even * (1.0f - ACO_SPACING_TOLERANCE)) return 1.0f;
    return std::clamp(gap / even, ACO_SPACING_MIN_SPEED, 1.0f);
}

// ----------------------------------------------------------------------------
// Where each ACO drone is along the tour
// ----------------------------------------------------------------------------
// A drone on the tour is flying from waypoint (target - 1) to waypoint target, so its arc position is that
// waypoint's arc position plus how far it has come from it (capped at the leg's length).
std::vector<std::pair<size_t, float>> MissionController::acoTourPositions() const {
    std::vector<std::pair<size_t, float>> out;
    const std::vector<Vec2>& tour = plan_.acoTour;
    const size_t n = tour.size();
    if (n < 2 || acoArc_.size() != n || acoTotal_ <= 0.0f) return out;

    for (size_t i = 0; i < drones_.size(); ++i) {
        const Drone& d = drones_[i];
        // Only drones actually flying the tour: not in transit to it, holding, returning or on another route
        if (!d.isACO || d.state != DroneState::Active || d.holding || !d.transit.empty() || d.route.size() != n) continue;
        const size_t prev = (d.target + n - 1) % n;
        const float leg = std::sqrt(sqDist(tour[prev], tour[d.target]));
        const float along = std::min(std::sqrt(sqDist(tour[prev], d.pos)), leg);
        out.emplace_back(i, std::fmod(acoArc_[prev] + along, acoTotal_));
    }
    std::sort(out.begin(), out.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
    return out;
}

// ----------------------------------------------------------------------------
// Gap from each ACO drone to the next one ahead (the tour is flown in increasing arc position)
// ----------------------------------------------------------------------------
std::vector<float> MissionController::acoGaps() const {
    const auto pos = acoTourPositions();
    std::vector<float> gaps;
    for (size_t k = 0; k < pos.size(); ++k) {
        const float ahead = pos[(k + 1) % pos.size()].second;
        float gap = ahead - pos[k].second;
        if (gap <= 0.0f) gap += acoTotal_;   // The drone ahead is past the tour's start (or it is the only drone)
        gaps.push_back(gap);
    }
    return gaps;
}

// ----------------------------------------------------------------------------
// Set this frame's speed of every ACO drone from its gap to the drone ahead
// ----------------------------------------------------------------------------
void MissionController::updateAcoSpacing() {
    for (Drone& d : drones_) d.speedScale = 1.0f;   // Everyone else (and anyone not on the tour) at full speed
    const auto pos = acoTourPositions();
    if (pos.size() < 2) return;
    const std::vector<float> gaps = acoGaps();
    const float even = acoTotal_ / static_cast<float>(pos.size());
    for (size_t k = 0; k < pos.size(); ++k) drones_[pos[k].first].speedScale = acoSpacingSpeed(gaps[k], even);
}
