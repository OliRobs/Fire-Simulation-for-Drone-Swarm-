#include "fire_model.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

// ============================================================================
// Wildfire Propagation: Rothermel rates on a cellular automaton
//
// 1. Per fuel model: the no-wind no-slope Rothermel fuel bed (reaction intensity, R_0, wind/slope
//    coefficients) for the current fuel moisture. Only 40 models exist, so this is done once per model.
// 2. Per cell: wind (20-ft wind x wind adjustment factor = midflame wind) and the cell's slope and aspect are
//    combined by vector addition into a head fire rate, direction and ellipse (RMRS-GTR-371 section 6).
// 3. Per cell and direction: the ellipse gives the rate of spread towards each of the 8 neighbours.
// 4. Fire arrival times spread outwards through the grid in time order (a priority queue), so the front
//    moves at exactly the Rothermel rates however long a tick is.
// ============================================================================

namespace {
    constexpr double kNever = std::numeric_limits<double>::infinity();
    constexpr int kCodeSlots = 256;   // FBFM40 numbers fit in 8 bits (RMRS-GTR-153 p. 5)
}

// ----------------------------------------------------------------------------
// Wind direction
// ----------------------------------------------------------------------------
// Meteorological convention: a "225 degree" wind comes FROM the south-west and blows towards the north-east.
// Grid frame: x = east, y = south (rows run down the map, north-up rasters), so north is -y.
Vec2 windTowardGrid(float fromDeg) {
    constexpr float DEG_TO_RAD = 0.0174532925f;
    const float toRad = (fromDeg + 180.0f) * DEG_TO_RAD;
    return Vec2(std::sin(toRad), -std::cos(toRad));
}

// ----------------------------------------------------------------------------
// Simulation set-up: precompute every cell's directional spread rates
// ----------------------------------------------------------------------------
FireSimulation::FireSimulation(std::vector<Cell> grid, const FireWeather& weather,
    float cellSizeX, float cellSizeY, float minutesPerTick)
    : grid_(std::move(grid)), fireMask_(grid_.size(), 0), minutesPerTick_(minutesPerTick) {
    constexpr int N = GRID_WIDTH * GRID_HEIGHT;

    // Unit direction (in metres, so non-square cells still give true directions) and length of each step
    double stepUx[NUM_DIRS], stepUy[NUM_DIRS];
    for (int d = 0; d < NUM_DIRS; ++d) {
        const double sx = DX[d] * cellSizeX, sy = DY[d] * cellSizeY;
        const double len = std::sqrt(sx * sx + sy * sy);
        stepDistance_[d] = static_cast<float>(len);
        stepUx[d] = sx / len;
        stepUy[d] = sy / len;
    }

    // --- 1. One Rothermel fuel bed per fuel model in use ---
    SurfaceFuelBed beds[kCodeSlots];
    double midflame[kCodeSlots] = {};          // Midflame wind for each model, ft/min
    bool known[kCodeSlots] = {};
    const double wind20 = weather.wind20ftKmh * KMH_TO_FT_PER_MIN;
    for (const Cell& c : grid_) {
        const int code = c.fuelModel;
        if (code <= 0 || code >= kCodeSlots || known[code]) continue;
        known[code] = true;
        if (const FuelModel* m = findFuelModel(code)) {
            beds[code] = computeFuelBed(*m, weather.moisture);
            // Fuel depth sets how much the vegetation slows the 20-ft wind, unless a WAF is given
            const double waf = weather.windAdjustmentFactor > 0.0f ? weather.windAdjustmentFactor
                                                                   : unshelteredWaf(m->depthFt);
            midflame[code] = wind20 * waf;
        }
    }
    const Vec2 wind = windTowardGrid(weather.windFromDeg);

    // --- 2 & 3. Per-cell ellipse and directional rates ---
    dirRos_.assign(static_cast<size_t>(N) * NUM_DIRS, 0.0f);
    headRos_.assign(N, 0.0f);
    burnDuration_.assign(N, 0.0f);
    burnable_.assign(N, 0);
    arrival_.assign(N, kNever);
    ignitionTime_.assign(N, kNever);

    constexpr double DEG_TO_RAD = 0.017453292519943295;
    const double minCell = std::min(cellSizeX, cellSizeY);
    for (int idx = 0; idx < N; ++idx) {
        const Cell& c = grid_[idx];
        const int code = c.fuelModel;
        if (code <= 0 || code >= kCodeSlots || !beds[code].burnable()) continue;   // No fuel, or too wet to burn
        if (ROADS_ARE_FUEL_BREAKS && c.roadClass >= ROAD_CLASS_MIN && c.roadClass <= ROAD_CLASS_MAX) continue;

        const SurfaceFuelBed& bed = beds[code];
        const SpreadEllipse e = computeSpread(bed, midflame[code], wind.x, wind.y,
            std::tan(c.slopeDeg * DEG_TO_RAD), c.upslopeX, c.upslopeY, weather.windLimit);

        burnable_[idx] = 1;
        headRos_[idx] = static_cast<float>(e.headRos * FT_PER_MIN_TO_M_PER_MIN);
        for (int d = 0; d < NUM_DIRS; ++d) {
            dirRos_[static_cast<size_t>(idx) * NUM_DIRS + d] =
                static_cast<float>(e.rosToward(stepUx[d], stepUy[d]) * FT_PER_MIN_TO_M_PER_MIN);
        }
        // A 30 m cell is "on fire" while the front crosses it plus the flaming residence time t_r (our
        // choice for the automaton: the Rothermel flame depth itself is much less than a cell).
        burnDuration_[idx] = static_cast<float>(minCell / headRos_[idx] + bed.residenceTime);
    }
}

// ----------------------------------------------------------------------------
// Ignite a cell at time t and schedule fire arrival at its neighbours
// ----------------------------------------------------------------------------
void FireSimulation::igniteAt(int idx, double t) {
    grid_[idx].state = CellState::Burning;
    fireMask_[idx] = 1;      // Binary fire mask for distance transform / spatial tracking
    dirty_.push_back(idx);   // Queue cell index for renderer updates
    fireGrew_ = true;        // Flag to trigger path planning or search graph invalidation
    ignitionTime_[idx] = t;
    arrival_[idx] = t;
    burning_.push_back(idx);

    const int x = idx % GRID_WIDTH, y = idx / GRID_WIDTH;
    for (int d = 0; d < NUM_DIRS; ++d) {
        const int nx = x + DX[d], ny = y + DY[d];
        if (nx < 0 || nx >= GRID_WIDTH || ny < 0 || ny >= GRID_HEIGHT) continue;
        const int nIdx = ny * GRID_WIDTH + nx;
        if (!burnable_[nIdx] || grid_[nIdx].state != CellState::Unburned) continue;

        // Diagonal steps may not squeeze between two nonburnable cells: otherwise a one-cell-wide diagonal
        // road or river would leak fire through its corners.
        if (DX[d] != 0 && DY[d] != 0 && !burnable_[y * GRID_WIDTH + nx] && !burnable_[ny * GRID_WIDTH + x]) continue;

        // Half the step at this cell's rate, half at the neighbour's: rate = harmonic mean of the two
        const float rOut = dirRos_[static_cast<size_t>(idx) * NUM_DIRS + d];
        const float rIn = dirRos_[static_cast<size_t>(nIdx) * NUM_DIRS + d];
        if (rOut <= 0.0f || rIn <= 0.0f) continue;
        const double travel = 0.5 * stepDistance_[d] / rOut + 0.5 * stepDistance_[d] / rIn;   // minutes

        const double tArrive = t + travel;
        if (tArrive < arrival_[nIdx]) {   // Keep only the earliest arrival; older queue entries go stale
            arrival_[nIdx] = tArrive;
            queue_.push({ tArrive, nIdx });
        }
    }
}

// ----------------------------------------------------------------------------
// Public External Ignition Trigger
// ----------------------------------------------------------------------------
// Manually ignites a cell (e.g., initial fire origin or lightning strikes) at the current time.
bool FireSimulation::ignite(int idx) {
    if (grid_[idx].state != CellState::Unburned) return false; // Already burning/burned
    if (!burnable_[idx]) return false;                         // Nonburnable (water, urban, road...) or too wet
    igniteAt(idx, now_);
    return true;
}

// ----------------------------------------------------------------------------
// Simulation Step (one tick of SIM_MINUTES_PER_TICK)
// ----------------------------------------------------------------------------
// 1. Advance the clock and ignite, in time order, every cell the fire reaches by then (a cell lit early in
//    the tick can light its own neighbours within the same tick if the fire is fast).
// 2. Burning cells whose burn duration has passed become Burned.
void FireSimulation::step() {
    now_ += minutesPerTick_;

    while (!queue_.empty() && queue_.top().first <= now_) {
        const Event ev = queue_.top();
        queue_.pop();
        const int idx = ev.second;
        // Skip stale entries: the cell already burned, or a faster route was found after this was queued
        if (grid_[idx].state != CellState::Unburned || ev.first > arrival_[idx]) continue;
        igniteAt(idx, ev.first);
    }

    size_t keep = 0; // In-place compaction of the burning list
    for (size_t i = 0; i < burning_.size(); ++i) {
        const int idx = burning_[i];
        Cell& c = grid_[idx];
        c.burnTime = static_cast<float>(now_ - ignitionTime_[idx]);
        if (c.burnTime >= burnDuration_[idx]) {
            c.state = CellState::Burned;
            dirty_.push_back(idx);
        }
        else {
            burning_[keep++] = idx;
        }
    }
    burning_.resize(keep);
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

// Checks and resets the fire growth flag
bool FireSimulation::takeFireGrew() {
    const bool grew = fireGrew_;
    fireGrew_ = false;
    return grew;
}
