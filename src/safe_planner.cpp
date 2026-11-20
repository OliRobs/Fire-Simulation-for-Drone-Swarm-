#include "safe_planner.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>
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

// Constructor: allocate the A* buffers once (one entry per grid cell)
SafePlanner::SafePlanner(const DistanceTransform& liveEdt)
    : edt_(liveEdt),
      nodes_(GRID_WIDTH * GRID_HEIGHT, Node{ 0.0f, -1, 0, 0 }) {
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
// Flyable-Cell Bitmap
// ----------------------------------------------------------------------------
// Brings fly_ up to date with the live distance field: fly_ at a cell is exactly cellOk() there, so A* finds the
// same paths whether it reads the bitmap or the field. Reading one byte (with no bounds checks, thanks to the
// border) is much cheaper than cellOk's four comparisons and an 8-byte load, and A* makes about 20 such tests
// per cell it expands. Rebuilt only when the field has changed since the last search.
void SafePlanner::refreshFlyable() {
    if (flyBuilt_ && edt_.version() == flyVersion_) return;   // Field unchanged: the bitmap is still right
    if (!flyBuilt_) {
        // First use: border cells 0, every grid cell flyable until the field says otherwise
        fly_.assign(static_cast<size_t>(FLY_W) * (GRID_HEIGHT + 2), 0);
        for (int y = 0; y < GRID_HEIGHT; ++y) {
            std::fill(fly_.begin() + (y + 1) * FLY_W + 1, fly_.begin() + (y + 1) * FLY_W + 1 + GRID_WIDTH, std::uint8_t(1));
        }
        flyBuilt_ = true;
        flyWin_ = GridRect{};
    }
    // Cells the previous window covered are back to flyable (outside the new window the field reads "very far")
    for (int y = flyWin_.y0; y < flyWin_.y1; ++y) {
        std::fill(fly_.begin() + (y + 1) * FLY_W + 1 + flyWin_.x0, fly_.begin() + (y + 1) * FLY_W + 1 + flyWin_.x1,
                  std::uint8_t(1));
    }
    // Inside the field's window: the same test as cellOk()
    const GridRect win = edt_.window();
    for (int y = win.y0; y < win.y1; ++y) {
        std::uint8_t* row = fly_.data() + (y + 1) * FLY_W + 1;
        for (int x = win.x0; x < win.x1; ++x) row[x] = edt_.squaredDistance(y * GRID_WIDTH + x) >= HARD_NO_FLY_SQ;
    }
    flyWin_ = win;
    flyVersion_ = edt_.version();
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
// Returns the direct line if it is clear; otherwise runs A* around the hard no-fly zone and
// applies string-pulling. Returns an EMPTY vector when no safe path exists (see header).
std::vector<Vec2> SafePlanner::leg(Vec2 a, Vec2 b) {
    // Fast path: direct line of sight is clear
    if (lineClear(a, b)) return { b };

    // Snap start and destination to the nearest flyable cells. If either has no flyable cell
    // nearby there is nothing safe to fly to, so report failure rather than a straight line.
    const int s = snap(a), t = snap(b);
    if (s < 0 || t < 0) return {};

    // Does b itself lie in flyable space? If not, the path will finish at the snapped cell instead.
    const bool goalFlyable = cellOk(static_cast<int>(std::floor(b.x)), static_cast<int>(std::floor(b.y)));

    // --- Phase 1: A* Grid Pathfinding ---
    refreshFlyable();   // fly_[...] below is exactly cellOk(), read from the bitmap
    const std::uint8_t* fly = fly_.data();
    // New search: bumping the stamp invalidates every node's g/came from previous searches at once.
    if (++search_ == 0) {   // Stamp wrapped round (after ~4 billion searches): reset properly
        for (Node& nd : nodes_) { nd.openStamp = 0; nd.closedStamp = 0; }
        search_ = 1;
    }
    Node* const nodes = nodes_.data();

    // Octile distance heuristic to the goal from cell (x, y), admissible for 8-connected moves:
    // |dx| + |dy| - (2 - sqrt 2) * min(|dx|, |dy|). Same float arithmetic as before, with the cell's coordinates
    // passed in rather than recovered from its index by a divide and a modulo.
    const int tx = t % GRID_WIDTH, ty = t / GRID_WIDTH;
    auto hTo = [tx, ty](int x, int y) {
        const float dx = static_cast<float>(std::abs(x - tx)), dy = static_cast<float>(std::abs(y - ty));
        return dx + dy - 0.58579f * std::min(dx, dy);
    };

    // --- Open list: buckets by f-cost (see safe_planner.h), each a min-heap of (f, cell) ---
    using QItem = std::pair<float, int>;   // {f_score, cell_index}
    const std::greater<QItem> after;       // Heap order: smallest (f, cell) at the front
    size_t lowest = 0;                     // No bucket below this one holds an entry
    size_t highest = 0;                    // Highest bucket used in this search (cleared at the end)
    size_t queued = 0;                     // Entries in all buckets
    auto push = [&](float f, int v) {
        const size_t k = static_cast<size_t>(f);   // f >= 0, so this is floor(f)
        if (k >= buckets_.size()) buckets_.resize(k + 1);
        buckets_[k].push_back({ f, v });
        std::push_heap(buckets_[k].begin(), buckets_[k].end(), after);
        lowest = std::min(lowest, k);              // (f never drops below the last pop, but be safe)
        highest = std::max(highest, k);
        ++queued;
    };
    auto pop = [&]() {
        while (buckets_[lowest].empty()) ++lowest;
        std::vector<QItem>& bk = buckets_[lowest];
        std::pop_heap(bk.begin(), bk.end(), after);
        const int u = bk.back().second;
        bk.pop_back();
        --queued;
        return u;
    };

    nodes[s].g = 0.0f;
    nodes[s].came = -1;
    nodes[s].openStamp = search_;
    {
        const float f0 = hTo(s % GRID_WIDTH, s / GRID_WIDTH);
        lowest = static_cast<size_t>(f0);
        push(f0, s);
    }

    bool reached = false;
    while (queued > 0) {
        const int u = pop();

        Node& nu = nodes[u];
        if (nu.closedStamp == search_) continue;  // Stale queue entry: already expanded
        nu.closedStamp = search_;

        if (u == t) { reached = true; break; }    // Goal reached

        const int ux = u % GRID_WIDTH, uy = u / GRID_WIDTH;
        const int fu = (uy + 1) * FLY_W + (ux + 1);   // u in the padded bitmap

        // Expand all 8 neighbours
        for (int d = 0; d < NUM_DIRS; ++d) {
            const int nx = ux + DX[d], ny = uy + DY[d];
            if (!fly[fu + DY[d] * FLY_W + DX[d]]) continue;   // cellOk(nx, ny); off the grid reads 0

            // No corner cutting: a diagonal step needs both orthogonal cells it squeezes past to be flyable
            const bool diagonal = DX[d] != 0 && DY[d] != 0;
            if (diagonal && (!fly[fu + DX[d]] || !fly[fu + DY[d] * FLY_W])) continue;

            const int v = ny * GRID_WIDTH + nx;
            const float ng = nu.g + (diagonal ? 1.41421f : 1.0f);

            // Better than the cost so far (infinity if v has not been reached in this search)?
            Node& nv = nodes[v];
            const float gv = nv.openStamp == search_ ? nv.g : std::numeric_limits<float>::max();
            if (ng < gv) {
                nv.g = ng;
                nv.came = u;
                nv.openStamp = search_;
                push(ng + hTo(nx, ny), v);
            }
        }
    }

    // Empty the buckets for the next search (entries are left when the goal is reached early)
    for (size_t k = 0; k <= highest && k < buckets_.size(); ++k) buckets_[k].clear();

    // The goal is cut off from the start by fire: there is no safe leg.
    if (!reached) return {};

    // --- Phase 2: Reconstruct the cell path (start point, cell centres..., end point) ---
    std::vector<Vec2> full;
    for (int v = t; v != -1; v = nodes[v].came) {
        full.push_back(cellCentre(v));
    }
    full.push_back(a);
    std::reverse(full.begin(), full.end());
    if (goalFlyable) full.push_back(b);   // Otherwise finish at the centre of the snapped cell

    // --- Phase 3: String pulling ---
    // From each point, jump to the farthest later point that is still in clear line of sight.
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
