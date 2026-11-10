#pragma once
#include <cstdint>
#include <functional>
#include <queue>
#include <utility>
#include <vector>
#include "crown_fire.h"
#include "rothermel.h"
#include "types.h"

// Weather inputs to the Rothermel model. None of these are in the LANDFIRE layers: the defaults are a
// fixed scenario (config.h, and Scott & Burgan's D2L2 moisture) until real weather data is wired in.
struct FireWeather {
    FuelMoisture moisture;                          // Dead 1/10/100-h and live herb/woody, fractions
    float wind20ftKmh = WIND_SPEED_20FT_KMH;        // Open wind at 20 ft
    float windFromDeg = WIND_FROM_DEG;              // Meteorological direction (wind comes FROM), clockwise from N
    float windAdjustmentFactor = WIND_ADJUSTMENT_FACTOR; // <= 0: from fuel depth and the canopy layers (Andrews 2012)
    bool windLimit = APPLY_WIND_LIMIT;
    bool crownFire = ENABLE_CROWN_FIRE;              // Allow torching / crown fire where there is canopy
    float foliarMoisturePct = FOLIAR_MOISTURE_PCT;   // Live canopy foliage moisture, %
};

// Unit vector, in the grid frame (x = east, y = south), that a wind FROM fromDeg blows towards
Vec2 windTowardGrid(float fromDeg);

// Rothermel-driven cellular automaton.
//
// Each cell gets a spread ellipse from its fuel model, slope, aspect and the wind (rothermel.h), with the wind
// reduced under a forest canopy. Where there is canopy and the surface fire is intense enough, the cell becomes
// a passive or active crown fire (crown_fire.h); active crown fire spreads at the crown rate. The rate at
// which fire crosses from a cell to a neighbour is the harmonic mean of the two cells' rates in that
// direction (half the step through each cell), so the travel time is exact for the Rothermel rates. Fire
// arrival times are then propagated with a priority queue (shortest travel time first), and step() advances
// the clock by SIM_MINUTES_PER_TICK, igniting every cell whose arrival time has come. This makes spread
// deterministic and independent of the tick length.
class FireSimulation {
public:
    // cellSizeX/Y: ground size of one cell in metres
    FireSimulation(std::vector<Cell> grid, const FireWeather& weather = FireWeather(),
                   float cellSizeX = DEFAULT_CELL_SIZE_M, float cellSizeY = DEFAULT_CELL_SIZE_M,
                   float minutesPerTick = SIM_MINUTES_PER_TICK);

    bool ignite(int idx);   // Manual ignition (start point, mouse click). False if the cell isn't Unburned or can't burn.
    void step();            // Advance one tick

    bool canBurn(int idx) const { return burnable_[idx] != 0; }  // Burnable fuel, not a road fuel break
    float headRos(int idx) const { return headRos_[idx]; }       // Head fire rate of spread, m/min
    double minutes() const { return now_; }                       // Simulated time since the start

    const std::vector<Cell>& grid() const { return grid_; }
    const std::vector<std::uint8_t>& fireMask() const { return fireMask_; } // 1 = Burning or Burned

    std::vector<int> takeDirtyCells(); // Cells whose state changed since the last call (for repainting)
    bool takeFireGrew();               // True if any cell ignited since the last call

private:
    void igniteAt(int idx, double t);  // Mark a cell alight at time t and schedule its neighbours

    std::vector<Cell> grid_;
    float stepDistance_[NUM_DIRS];     // Metres between cell centres in each direction
    std::vector<float> dirRos_;        // [cell * 8 + direction]: Rothermel rate from the cell in that direction, m/min
    std::vector<float> headRos_;       // Head fire rate per cell, m/min (0 = can't burn)
    std::vector<float> burnDuration_;  // Minutes a cell shows as Burning: front crossing time + flame residence time
    std::vector<std::uint8_t> burnable_;
    std::vector<double> arrival_;      // Earliest scheduled fire arrival time per cell (minutes; +inf = none)
    std::vector<double> ignitionTime_;

    using Event = std::pair<double, int>;  // (arrival time, cell)
    std::priority_queue<Event, std::vector<Event>, std::greater<Event>> queue_;

    std::vector<int> burning_;         // Cells currently Burning
    std::vector<std::uint8_t> fireMask_;
    std::vector<int> dirty_;
    bool fireGrew_ = false;
    double now_ = 0.0;
    double minutesPerTick_;
};
