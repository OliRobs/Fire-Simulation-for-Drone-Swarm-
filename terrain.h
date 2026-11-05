#pragma once
#include <cstdint>
#include <random>
#include <vector>
#include "types.h"

// Builds the grid from the DEM and fuel GeoTIFFs, with procedural fallbacks if a file is missing.
std::vector<Cell> loadTerrain(const char* demPath, const char* fuelPath, std::mt19937& rng);

// Elevation mapped linearly to 0 (lowest, black) .. 255 (highest, white), one value per cell.
std::vector<std::uint8_t> computeTerrainShade(const std::vector<Cell>& grid);
