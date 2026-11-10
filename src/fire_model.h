#pragma once
#include <cstdint>
#include <random>
#include <vector>
#include "types.h"

// --- Physics ---
float getWindFactor(int dx, int dy, const Vec2& wind);
// Multiplier for spreading from a cell to a neighbour `horizontalDistance` away (same units as elevation).
float getSlopeFactor(float elevationCurrent, float elevationNeighbor, float horizontalDistance);

// Cellular-automaton wildfire. Only burning cells are visited each tick, and ignition probabilities
// (wind, slope and the neighbour's fuel are static until a neighbour burns) are precomputed once.
class FireSimulation {
public:
    // cellSizeX/Y: ground size of one cell in metres, so slope is rise over real distance
    FireSimulation(std::vector<Cell> grid, Vec2 wind, std::uint32_t seed,
                   float cellSizeX = DEFAULT_CELL_SIZE_M, float cellSizeY = DEFAULT_CELL_SIZE_M);

    bool ignite(int idx);   // Manual ignition (start point, mouse click). False if the cell isn't Unburned or can't burn.
    void step();            // Advance one tick

    const std::vector<Cell>& grid() const { return grid_; }
    const std::vector<std::uint8_t>& fireMask() const { return fireMask_; } // 1 = Burning or Burned

    std::vector<int> takeDirtyCells(); // Cells whose state changed since the last call (for repainting)
    bool takeFireGrew();               // True if any cell ignited since the last call

private:
    void igniteCell(int idx);

    std::vector<Cell> grid_;
    std::vector<float> ignitionProb_;  // [cell * 8 + direction]; 0 for out-of-bounds or fuel-less neighbours
    std::vector<int> burning_;         // Active frontier
    std::vector<int> newlyIgnited_;
    std::vector<std::uint8_t> fireMask_;
    std::vector<int> dirty_;
    bool fireGrew_ = false;
    std::mt19937 rng_;
    std::uniform_real_distribution<float> dist_{ 0.0f, 1.0f };
};
