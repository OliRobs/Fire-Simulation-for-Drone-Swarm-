#pragma once
#include <cstdint>
#include <random>
#include <vector>
#include "types.h"

// The simulation grid plus the ground size of one cell, which the slope calculation needs.
struct Terrain {
    std::vector<Cell> grid;
    float cellSizeX = DEFAULT_CELL_SIZE_M;   // Metres per cell, east-west
    float cellSizeY = DEFAULT_CELL_SIZE_M;   // Metres per cell, north-south
};

// Builds the grid from the DEM and fuel GeoTIFFs, with procedural fallbacks if a file is missing.
// The DEM defines the grid (coordinate system and extent); the fuel map is warped onto it. Fuel maps holding
// FBFM40 codes are converted through the Scott & Burgan lookup; maps already in 0..1 are used as flammability.
Terrain loadTerrain(const char* demPath, const char* fuelPath, std::mt19937& rng);

// Elevation mapped linearly to 0 (lowest, black) .. 255 (highest, white), one value per cell.
std::vector<std::uint8_t> computeTerrainShade(const std::vector<Cell>& grid);
