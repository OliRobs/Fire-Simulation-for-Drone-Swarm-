#include "fire_model.h"
#include <algorithm>
#include <cmath>
#include <utility>

// ============================================================================
// Wildfire Propagation Physics & Simulation Model
//
// Models 2D grid-based fire spread influenced by:
// 1. Wind speed and direction (vector dot product modifier).
// 2. Terrain slope gradient (uphill elevation delta exponential multiplier).
// 3. Cell fuel load and ignition probability per neighbor direction.
// ============================================================================

// ----------------------------------------------------------------------------
// Wind Influence Calculation
// ----------------------------------------------------------------------------
// Returns an exponential multiplier based on how closely aligned the spread 
// direction (dx, dy) is with the ambient wind direction vector.
float getWindFactor(int dx, int dy, const Vec2& wind) {
    const float len = std::sqrt(static_cast<float>(dx * dx + dy * dy));
    if (len == 0.0f) return 1.0f; // Center cell (no directional step)

    // Dot product between normalized direction and wind vector
    const float dot = (dx / len) * wind.x + (dy / len) * wind.y;

    // Exponential scaling: downwind spreads significantly faster
    return std::exp(dot * 0.8f);
}

// ----------------------------------------------------------------------------
// Slope / Elevation Influence Calculation
// ----------------------------------------------------------------------------
// Returns an exponential multiplier based on elevation differences.
// Fires spread exponentially faster when traveling uphill to neighboring cells.
float getSlopeFactor(float elevationCurrent, float elevationNeighbor) {
    return std::exp((elevationNeighbor - elevationCurrent) * 0.05f);
}

// ----------------------------------------------------------------------------
// Simulation Engine Initialization
// ----------------------------------------------------------------------------
// Precomputes static ignition probability lookup tables for all cells across 
// all 8 (or NUM_DIRS) neighbor direction vectors.
FireSimulation::FireSimulation(std::vector<Cell> grid, Vec2 wind, std::uint32_t seed)
    : grid_(std::move(grid)), fireMask_(grid_.size(), 0), rng_(seed) {
    constexpr int N = GRID_WIDTH * GRID_HEIGHT;

    // Precalculate directional wind factors for all neighbor directions
    float windFactor[NUM_DIRS];
    for (int d = 0; d < NUM_DIRS; ++d) {
        windFactor[d] = getWindFactor(DX[d], DY[d], wind);
    }

    // Allocate directional ignition probability array: N cells * NUM_DIRS directions
    ignitionProb_.assign(static_cast<size_t>(N) * NUM_DIRS, 0.0f);

    for (int y = 0; y < GRID_HEIGHT; ++y) {
        for (int x = 0; x < GRID_WIDTH; ++x) {
            const int idx = y * GRID_WIDTH + x;

            for (int d = 0; d < NUM_DIRS; ++d) {
                const int nx = x + DX[d], ny = y + DY[d];
                // Skip out-of-bounds neighbors
                if (nx < 0 || nx >= GRID_WIDTH || ny < 0 || ny >= GRID_HEIGHT) continue;

                const Cell& nb = grid_[ny * GRID_WIDTH + nx];
                if (nb.fuel <= 0.0f) continue; // Unburnable neighbor cell

                // Combined probability formula: BaseRate * WindMod * SlopeMod * FuelDensity
                const float p = BASE_IGNITION_RATE * windFactor[d]
                    * getSlopeFactor(grid_[idx].elevation, nb.elevation) * nb.fuel;

                // Store directional probability clamped to max 1.0
                ignitionProb_[static_cast<size_t>(idx) * NUM_DIRS + d] = std::min(p, 1.0f);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Internal Cell State Ignition
// ----------------------------------------------------------------------------
// Updates internal state when a cell transition to Burning occurs.
void FireSimulation::igniteCell(int idx) {
    grid_[idx].state = CellState::Burning;
    fireMask_[idx] = 1;      // Binary fire mask for distance transform / spatial tracking
    dirty_.push_back(idx);   // Queue cell index for renderer/GPU updates
    fireGrew_ = true;        // Flag to trigger path planning or search graph invalidation
}

// ----------------------------------------------------------------------------
// Public External Ignition Trigger
// ----------------------------------------------------------------------------
// Manually ignites a cell (e.g., initial fire origin or lightning strikes).
bool FireSimulation::ignite(int idx) {
    if (grid_[idx].state != CellState::Unburned) return false; // Already burning/burned
    igniteCell(idx);
    burning_.push_back(idx); // Add cell to active burning frontier list
    return true;
}

// ----------------------------------------------------------------------------
// Simulation Step (Single Time Frame / Tick)
// ----------------------------------------------------------------------------
// Advances burning processes, calculates stochastic neighbor ignitions, 
// consumes fuel, and transitions fully consumed cells to Burned state.
void FireSimulation::step() {
    newlyIgnited_.clear();
    size_t keep = 0; // In-place compaction index for active burning frontier

    // Process all currently burning frontier cells
    for (size_t i = 0; i < burning_.size(); ++i) {
        const int idx = burning_[i];
        const float* probs = &ignitionProb_[static_cast<size_t>(idx) * NUM_DIRS];

        // Attempt stochastic ignition of adjacent unburned cells
        for (int d = 0; d < NUM_DIRS; ++d) {
            const float p = probs[d];
            if (p <= 0.0f) continue;

            const int nIdx = idx + DY[d] * GRID_WIDTH + DX[d];
            if (grid_[nIdx].state != CellState::Unburned) continue;

            // Monte Carlo trial based on precomputed probability
            if (dist_(rng_) < p) {
                igniteCell(nIdx);
                newlyIgnited_.push_back(nIdx);
            }
        }

        // Consume fuel of current burning cell
        Cell& c = grid_[idx];
        c.fuel -= FUEL_BURN_RATE;
        c.burnTime += 1.0f;

        if (c.fuel <= 0.0f) {
            // Fuel depleted: transition cell state to Burned
            c.state = CellState::Burned;
            dirty_.push_back(idx);
        }
        else {
            // Cell still has fuel: retain on burning frontier list
            burning_[keep++] = idx;
        }
    }

    // Shrink active list to remaining burning cells and append newly ignited cells
    burning_.resize(keep);
    burning_.insert(burning_.end(), newlyIgnited_.begin(), newlyIgnited_.end());
}

// ----------------------------------------------------------------------------
// State Tracking Queries
// ----------------------------------------------------------------------------

// Flushes and returns grid indices modified since the last call (for rendering optimization)
std::vector<int> FireSimulation::takeDirtyCells() {
    std::vector<int> out;
    out.swap(dirty_);
    return out;
}

// Atomically checks and resets the fire perimeter expansion flag
bool FireSimulation::takeFireGrew() {
    const bool grew = fireGrew_;
    fireGrew_ = false;
    return grew;
}