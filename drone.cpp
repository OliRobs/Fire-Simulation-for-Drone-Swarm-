#include "drone.h"
#include <algorithm>
#include <cmath>
#include <limits>

// ============================================================================
// Drone Navigation and Safety System
// 
// Manages path following, transit routing, and emergency fire-avoidance physics
// for autonomous aerial coverage/monitoring drones.
// ============================================================================

// ----------------------------------------------------------------------------
// Find Nearest Waypoint Index
// ----------------------------------------------------------------------------
// Searches a 2D route trajectory to find the index of the waypoint closest to 
// a given point `p`. Used for route matching or initial target selection.
size_t nearestWaypoint(const std::vector<Vec2>& route, Vec2 p) {
    size_t best = 0;
    float bd = std::numeric_limits<float>::max();

    for (size_t i = 0; i < route.size(); ++i) {
        const float dd = sqDist(route[i], p);
        if (dd < bd) {
            bd = dd;
            best = i;
        }
    }
    return best;
}

// ----------------------------------------------------------------------------
// Advance Drone Position Along Route/Transit (Fixed Speed Step)
// ----------------------------------------------------------------------------
// Moves the drone forward at `DRONE_SPEED` per tick. 
// Priority order:
// 1. Clears any queued fire-safe transit waypoints first.
// 2. Cycles continuously around its assigned circular coverage route.
void stepDrone(Drone& d) {
    float remaining = DRONE_SPEED; // Remaining movement distance budget for this frame
    const size_t maxSteps = d.route.size() + d.transit.size() + 1; // Safety guard against infinite loops

    for (size_t guard = 0; remaining > 1e-4f && guard <= maxSteps; ++guard) {
        const bool inTransit = !d.transit.empty();

        // Return if drone has finished transit and has no loop route assigned
        if (!inTransit && d.route.empty()) return;

        // Current target is the front of the transit queue, or the next point in the loop route
        const Vec2 goal = inTransit ? d.transit.front() : d.route[d.target];
        const Vec2 delta = goal - d.pos;
        const float len = std::sqrt(delta.x * delta.x + delta.y * delta.y);

        if (len <= remaining) {
            // Drone reaches current target waypoint this step
            d.pos = goal;
            remaining -= len; // Deduct used distance budget

            if (inTransit) {
                // Pop reached point from temporary transit queue
                d.transit.erase(d.transit.begin());
            }
            else {
                // Advance to next target waypoint in circular loop
                d.target = (d.target + 1) % d.route.size();
                if (d.route.size() == 1) break; // Single waypoint: hover in place
            }
        }
        else {
            // Move partial distance toward current waypoint and consume remaining budget
            d.pos += delta * (remaining / len);
            remaining = 0.0f;
        }
    }
}

// ----------------------------------------------------------------------------
// Emergency Fire Repulsion Safety Net
// ----------------------------------------------------------------------------
// Dynamic safety fallback: Fire expands dynamically between full path re-plans.
// Evaluates the Euclidean Distance Transform (EDT) gradient around the drone's position.
// If the drone is inside `HARD_MIN_FIRE_DISTANCE`, it is actively pushed outward 
// along the distance-field gradient away from the fire front.
void keepClearOfFire(Drone& d, const DistanceTransform& liveEdt) {
    const float minDist = static_cast<float>(std::sqrt(HARD_NO_FLY_SQ));

    // Helper lambda to look up Euclidean distance to the nearest fire front at grid (x, y)
    auto distAt = [&](int x, int y) {
        x = std::clamp(x, 0, GRID_WIDTH - 1);
        y = std::clamp(y, 0, GRID_HEIGHT - 1);
        return static_cast<float>(std::sqrt(liveEdt.squaredDistance(y * GRID_WIDTH + x)));
        };

    // Iterative gradient push-back solver (up to 40 step refinement passes)
    for (int iter = 0; iter < 40; ++iter) {
        const int x = std::clamp(static_cast<int>(d.pos.x), 0, GRID_WIDTH - 1);
        const int y = std::clamp(static_cast<int>(d.pos.y), 0, GRID_HEIGHT - 1);
        const float cur = distAt(x, y);

        // Safe distance satisfied; no intervention needed
        if (cur >= minDist) return;

        // Compute central difference spatial gradient of distance field vector: (dD/dx, dD/dy)
        Vec2 g(distAt(x + 1, y) - distAt(x - 1, y), distAt(x, y + 1) - distAt(x, y - 1));
        float gl = std::sqrt(g.x * g.x + g.y * g.y);

        // Handle zero-gradient flat regions (deep inside fire envelope)
        if (gl < 1e-3f) {
            // Fallback direction: Push drone radially outward from map center
            g = d.pos - Vec2(GRID_WIDTH * 0.5f, GRID_HEIGHT * 0.5f);
            gl = std::sqrt(g.x * g.x + g.y * g.y);

            // Extreme edge-case fallback: Default to positive X displacement
            if (gl < 1e-3f) { g = Vec2(1.0f, 0.0f); gl = 1.0f; }
        }

        // Repel drone along normalized gradient direction out to safe buffer distance (+0.3f padding)
        d.pos += g * ((minDist - cur + 0.3f) / gl);

        // Clamp final position within operational map grid boundaries
        d.pos.x = std::clamp(d.pos.x, 0.5f, GRID_WIDTH - 0.5f);
        d.pos.y = std::clamp(d.pos.y, 0.5f, GRID_HEIGHT - 0.5f);
    }
}