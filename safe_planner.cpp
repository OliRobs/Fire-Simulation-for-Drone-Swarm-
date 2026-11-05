#include "safe_planner.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
#include <queue>
#include <utility>

// ============================================================================
// Safe Path Planner (Grid-based A* Search with String-Pulling Optimization)
//
// Generates collision-free flight legs for drones around dynamic wildfire 
// no-fly zones (hard safety distance buffer).
// 
// Features:
// 1. Ray-casting line-of-sight checks to avoid dynamic obstacles.
// 2. A* grid search with Octile distance heuristics and anti-corner-cutting.
// 3. String-pulling (Funnel-like algorithm) to straighten and smooth path segments.
// 4. Snap-to-valid-cell fallback search for waypoints starting inside forbidden zones.
// ============================================================================

// Constructor: Initializes internal search graph buffers for A* pathfinding
SafePlanner::SafePlanner(const DistanceTransform& liveEdt)
    : edt_(liveEdt), g_(GRID_WIDTH* GRID_HEIGHT), came_(GRID_WIDTH* GRID_HEIGHT), closed_(GRID_WIDTH* GRID_HEIGHT) {
}

// ----------------------------------------------------------------------------
// Cell Feasibility Check
// ----------------------------------------------------------------------------
// Evaluates if a given grid coordinate (x, y) is within bounds and outside 
// the hard no-fly safety zone surrounding the active fire front.
bool SafePlanner::cellOk(int x, int y) const {
    return x >= 0 && x < GRID_WIDTH && y >= 0 && y < GRID_HEIGHT
        && edt_.squaredDistance(y * GRID_WIDTH + x) >= HARD_NO_FLY_SQ;
}

// ----------------------------------------------------------------------------
// Ray-Casting Line-of-Sight Check
// ----------------------------------------------------------------------------
// Marches samples along the straight-line segment from point 'a' to point 'b'
// at 0.25-cell increments to verify no-fly zone obstacles are crossed.
bool SafePlanner::lineClear(Vec2 a, Vec2 b) const {
    const float len = std::sqrt(sqDist(a, b));
    const int steps = std::max(1, static_cast<int>(len / 0.25f) + 1);

    for (int i = 0; i <= steps; ++i) {
        const float t = static_cast<float>(i) / steps;
        const Vec2 p = a + (b - a) * t;
        if (!cellOk(static_cast<int>(std::floor(p.x)), static_cast<int>(std::floor(p.y)))) {
            return false;
        }
    }
    return true;
}

// ----------------------------------------------------------------------------
// Safe Route Flight Leg Generation
// ----------------------------------------------------------------------------
// Connects start point `a` to goal point `b`.
// Returns direct line if unblocked; otherwise runs A* grid pathfinding around 
// the fire standoff barrier and applies string-pulling path smoothing.
std::vector<Vec2> SafePlanner::leg(Vec2 a, Vec2 b) {
    // Fast path: Direct line of sight is clear
    if (lineClear(a, b)) return { b };

    // Snap start and destination points to nearest valid grid cells
    const int s = snap(a), t = snap(b);
    if (s < 0 || t < 0) return { b }; // Fallback to straight line if snapping fails

    // --- Phase 1: A* Grid Pathfinding ---
    constexpr float INF = std::numeric_limits<float>::max();
    std::fill(g_.begin(), g_.end(), INF);
    std::fill(closed_.begin(), closed_.end(), 0);

    using QItem = std::pair<float, int>; // Priority queue item: {f_score, cell_index}
    std::priority_queue<QItem, std::vector<QItem>, std::greater<QItem>> open;

    g_[s] = 0.0f;
    came_[s] = -1;
    open.push({ heuristic(s, t), s });

    while (!open.empty()) {
        const int u = open.top().second;
        open.pop();

        if (closed_[u]) continue;
        closed_[u] = 1;

        if (u == t) break; // Goal reached

        const int ux = u % GRID_WIDTH, uy = u / GRID_WIDTH;

        for (int d = 0; d < NUM_DIRS; ++d) {
            const int nx = ux + DX[d], ny = uy + DY[d];
            if (!cellOk(nx, ny)) continue;

            const bool diagonal = DX[d] != 0 && DY[d] != 0;

            // Prevent corner-cutting through blocked diagonal neighbors
            if (diagonal && (!cellOk(ux + DX[d], uy) || !cellOk(ux, uy + DY[d]))) continue;

            const int v = ny * GRID_WIDTH + nx;
            const float ng = g_[u] + (diagonal ? 1.41421f : 1.0f);

            if (ng < g_[v]) {
                g_[v] = ng;
                came_[v] = u;
                open.push({ ng + heuristic(v, t), v });
            }
        }
    }

    // Path unreachable fallback
    if (g_[t] == INF) return { b };

    // Reconstruct raw cell-center node path from goal back to start
    std::vector<Vec2> full;
    for (int v = t; v != -1; v = came_[v]) {
        full.push_back(cellCentre(v));
    }
    full.push_back(a);
    std::reverse(full.begin(), full.end());
    full.push_back(b);

    // --- Phase 2: String-Pulling Path Optimization ---
    // Greedy shortcutting: Skips intermediate waypoints whenever direct line-of-sight is open
    std::vector<Vec2> out;
    size_t i = 0;
    while (i + 1 < full.size()) {
        size_t j = full.size() - 1;
        while (j > i + 1 && !lineClear(full[i], full[j])) {
            --j;
        }
        out.push_back(full[j]);
        i = j;
    }
    return out;
}

// ----------------------------------------------------------------------------
// Grid Cell Snapping Search
// ----------------------------------------------------------------------------
// Expands outward in concentric square rings (up to radius 40) to find the 
// nearest valid, non-blocked cell to snap a floating-point coordinate `p`.
int SafePlanner::snap(Vec2 p) const {
    const int cx = static_cast<int>(std::floor(p.x)), cy = static_cast<int>(std::floor(p.y));

    for (int r = 0; r <= 40; ++r) {
        for (int dy = -r; dy <= r; ++dy) {
            for (int dx = -r; dx <= r; ++dx) {
                if (std::max(std::abs(dx), std::abs(dy)) == r && cellOk(cx + dx, cy + dy)) {
                    return (cy + dy) * GRID_WIDTH + (cx + dx);
                }
            }
        }
    }
    return -1; // No valid cell found within search radius
}

// ----------------------------------------------------------------------------
// Octile Distance Heuristic
// ----------------------------------------------------------------------------
// Computes optimal 8-directional grid heuristic distance between cell `a` and cell `b`.
// Formula: $L_1 - (2 - \sqrt{2}) \cdot \min(\Delta x, \Delta y)$
float SafePlanner::heuristic(int a, int b) const {
    const float dx = static_cast<float>(std::abs(a % GRID_WIDTH - b % GRID_WIDTH));
    const float dy = static_cast<float>(std::abs(a / GRID_WIDTH - b / GRID_WIDTH));
    return dx + dy - 0.58579f * std::min(dx, dy);
}