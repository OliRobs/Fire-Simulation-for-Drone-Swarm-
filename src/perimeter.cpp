#include "perimeter.h"
#include "viewport.h"   // maskBounds

// ============================================================================
// Fire Perimeter Tracking & Containment System
//
// Manages the "committed" containment perimeter around the expanding wildfire.
// Instead of advancing continuously with every minor fire spark, the perimeter 
// updates discretely in local "sections" when the live fire encroaches within 
// a specific trigger threshold distance.
// ============================================================================

// Constructor: Initializes committed fire mask and computes the baseline 
// Euclidean Distance Transform (EDT) around the initial fire shape.
Perimeter::Perimeter(const std::vector<std::uint8_t>& initialFire)
    : committed_(initialFire), trigger_(initialFire.size(), 0),
    // Range limits (distance_transform.h): the committed-fire distance is read out to SEARCH_DISTANCE + 0.5
    // (mission zones; the ring and the trigger distance are nearer), the trigger distance only out to
    // PERIMETER_SECTION_RADIUS, so neither transform needs the whole grid.
    perimEdt_(GRID_WIDTH, GRID_HEIGHT, SEARCH_DISTANCE + 1),
    triggerEdt_(GRID_WIDTH, GRID_HEIGHT, PERIMETER_SECTION_RADIUS + 1) {
    // Initial distance field computation around committed fire boundaries
    perimEdt_.compute(committed_);
    // Record where the starting fire is, so later updates know where committed_ can be non-zero
    int x0, y0, x1, y1;
    if (maskBounds(committed_, GRID_WIDTH, GRID_HEIGHT, x0, y0, x1, y1)) fireSeen_ = GridRect{ x0, y0, x1 + 1, y1 + 1 };
}

// ----------------------------------------------------------------------------
// Update without a known fire rectangle: find it, then do the normal update
// ----------------------------------------------------------------------------
bool Perimeter::update(const std::vector<std::uint8_t>& liveFire) {
    int x0, y0, x1, y1;
    if (!maskBounds(liveFire, GRID_WIDTH, GRID_HEIGHT, x0, y0, x1, y1)) return update(liveFire, GridRect{});
    return update(liveFire, GridRect{ x0, y0, x1 + 1, y1 + 1 });
}

// ----------------------------------------------------------------------------
// Update Containment Perimeter
// ----------------------------------------------------------------------------
// Evaluates live fire expansion against the committed perimeter distance field.
// Returns `true` if any section of the committed perimeter was pushed forward.
bool Perimeter::update(const std::vector<std::uint8_t>& liveFire, const GridRect& fireBounds) {
    // Trigger Threshold: A live fire cell triggers a perimeter update once it 
    // encroaches within PERIMETER_UPDATE_DISTANCE of the current perimeter ring.
    constexpr double triggerSq = static_cast<double>(PERIMETER_DISTANCE - PERIMETER_UPDATE_DISTANCE)
        * (PERIMETER_DISTANCE - PERIMETER_UPDATE_DISTANCE);

    // Section Expansion Radius: Maximum influence distance around triggering cells
    constexpr double sectionSq = static_cast<double>(PERIMETER_SECTION_RADIUS) * PERIMETER_SECTION_RADIUS;

    // Every cell that is (or ever was) live fire, committed or a trigger lies inside 'area', so the loops
    // below only visit it: outside it liveFire, committed_ and trigger_ are all 0. Using the union with the
    // earlier rectangles (not just this one) also clears trigger_ cells set by a previous call.
    fireSeen_ = fireSeen_.unite(fireBounds);
    const GridRect area = fireSeen_;

    // --- Phase 1: Identify all live fire cells breaching the trigger distance ---
    bool any = false;
    for (int y = area.y0; y < area.y1; ++y) {
        for (int idx = y * GRID_WIDTH + area.x0, end = y * GRID_WIDTH + area.x1; idx < end; ++idx) {
            const bool trig = liveFire[idx] && perimEdt_.squaredDistance(idx) >= triggerSq;
            trigger_[idx] = trig;
            any = any || trig;
        }
    }

    // No cells breached the trigger distance; perimeter remains unchanged
    if (!any) return false;

    // --- Phase 2: Pull the committed perimeter forward around triggered sections ---
    // Compute distance field from the triggering cells (all inside 'area', so no grid scan is needed)
    triggerEdt_.compute(trigger_, area);

    for (int y = area.y0; y < area.y1; ++y) {
        for (int idx = y * GRID_WIDTH + area.x0, end = y * GRID_WIDTH + area.x1; idx < end; ++idx) {
            // Expand committed fire mask to include live fire cells within the section radius
            if (liveFire[idx] && !committed_[idx] && triggerEdt_.squaredDistance(idx) <= sectionSq) {
                committed_[idx] = 1;
            }
        }
    }

    // Recompute committed distance field following perimeter updates (committed_ lies inside 'area')
    perimEdt_.compute(committed_, area);
    return true;
}

// ----------------------------------------------------------------------------
// Extract Perimeter Ring Cells
// ----------------------------------------------------------------------------
// Returns a list of grid cell indices forming the visual/functional perimeter 
// ring band located at standard `PERIMETER_DISTANCE` away from committed fire.
std::vector<int> Perimeter::ringCells() const {
    constexpr double ringLo = (PERIMETER_DISTANCE - 0.5) * (PERIMETER_DISTANCE - 0.5);
    constexpr double ringHi = (PERIMETER_DISTANCE + 0.5) * (PERIMETER_DISTANCE + 0.5);

    // The ring is nearer than the transform's range (SEARCH_DISTANCE), so it lies inside the transform's
    // window: only that needs scanning, in the same row-major order as a whole-grid scan
    const GridRect win = perimEdt_.window();
    std::vector<int> ring;
    for (int y = win.y0; y < win.y1; ++y) {
        for (int idx = y * GRID_WIDTH + win.x0, end = y * GRID_WIDTH + win.x1; idx < end; ++idx) {
            const double d2 = perimEdt_.squaredDistance(idx);
            // Collect cells residing within the thin ring shell radius [PERIMETER_DISTANCE ± 0.5]
            if (d2 >= ringLo && d2 < ringHi) {
                ring.push_back(idx);
            }
        }
    }
    return ring;
}