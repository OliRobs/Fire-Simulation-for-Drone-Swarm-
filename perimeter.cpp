#include "perimeter.h"

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
    perimEdt_(GRID_WIDTH, GRID_HEIGHT), triggerEdt_(GRID_WIDTH, GRID_HEIGHT) {
    // Initial distance field computation around committed fire boundaries
    perimEdt_.compute(committed_);
}

// ----------------------------------------------------------------------------
// Update Containment Perimeter
// ----------------------------------------------------------------------------
// Evaluates live fire expansion against the committed perimeter distance field.
// Returns `true` if any section of the committed perimeter was pushed forward.
bool Perimeter::update(const std::vector<std::uint8_t>& liveFire) {
    constexpr int N = GRID_WIDTH * GRID_HEIGHT;

    // Trigger Threshold: A live fire cell triggers a perimeter update once it 
    // encroaches within PERIMETER_UPDATE_DISTANCE of the current perimeter ring.
    constexpr double triggerSq = static_cast<double>(PERIMETER_DISTANCE - PERIMETER_UPDATE_DISTANCE)
        * (PERIMETER_DISTANCE - PERIMETER_UPDATE_DISTANCE);

    // Section Expansion Radius: Maximum influence distance around triggering cells
    constexpr double sectionSq = static_cast<double>(PERIMETER_SECTION_RADIUS) * PERIMETER_SECTION_RADIUS;

    // --- Phase 1: Identify all live fire cells breaching the trigger distance ---
    bool any = false;
    for (int idx = 0; idx < N; ++idx) {
        const bool trig = liveFire[idx] && perimEdt_.squaredDistance(idx) >= triggerSq;
        trigger_[idx] = trig;
        any = any || trig;
    }

    // No cells breached the trigger distance; perimeter remains unchanged
    if (!any) return false;

    // --- Phase 2: Pull the committed perimeter forward around triggered sections ---
    // Compute distance field from the triggering cells
    triggerEdt_.compute(trigger_);

    for (int idx = 0; idx < N; ++idx) {
        // Expand committed fire mask to include live fire cells within the section radius
        if (liveFire[idx] && !committed_[idx] && triggerEdt_.squaredDistance(idx) <= sectionSq) {
            committed_[idx] = 1;
        }
    }

    // Recompute committed distance field following perimeter updates
    perimEdt_.compute(committed_);
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

    std::vector<int> ring;
    for (int idx = 0; idx < GRID_WIDTH * GRID_HEIGHT; ++idx) {
        const double d2 = perimEdt_.squaredDistance(idx);
        // Collect cells residing within the thin ring shell radius [PERIMETER_DISTANCE ± 0.5]
        if (d2 >= ringLo && d2 < ringHi) {
            ring.push_back(idx);
        }
    }
    return ring;
}