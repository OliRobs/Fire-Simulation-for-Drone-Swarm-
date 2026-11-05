#include "terrain.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include "raster_io.h"

// ============================================================================
// Terrain and Elevation/Fuel Data Manager
//
// Responsible for loading GIS raster datasets (GeoTIFFs) into the 2D grid,
// applying fallback synthetic terrain generation if files are missing, and
// calculating terrain hillshading for visualization.
// ============================================================================

// ----------------------------------------------------------------------------
// Load Terrain Data
// ----------------------------------------------------------------------------
// Attempts to load real-world elevation (DEM) and fuel map rasters via GDAL.
// If loading fails, it falls back to a synthetic procedural gradient terrain and 
// randomized fuel load.
std::vector<Cell> loadTerrain(const char* demPath, const char* fuelPath, std::mt19937& rng) {
    constexpr int N = GRID_WIDTH * GRID_HEIGHT;
    std::vector<float> rawElevation, rawFuel;

    std::cout << "Loading real-world GIS datasets via GDAL...\n";

    // Read elevation (demPath) and fuel (fuelPath) rasters resampled to grid dimensions
    const bool demLoaded = loadGeoTIFF(demPath, rawElevation, GRID_WIDTH, GRID_HEIGHT, true);
    const bool fuelLoaded = loadGeoTIFF(fuelPath, rawFuel, GRID_WIDTH, GRID_HEIGHT, false);

    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    std::vector<Cell> grid(N);

    for (int y = 0; y < GRID_HEIGHT; ++y) {
        for (int x = 0; x < GRID_WIDTH; ++x) {
            const int idx = y * GRID_WIDTH + x;

            // Assign elevation: use GeoTIFF data if loaded, else fallback to a synthetic diagonal slope
            grid[idx].elevation = demLoaded ? rawElevation[idx]
                : (x * 0.3f) + ((GRID_HEIGHT - y) * 0.5f);

            // Assign fuel density: use GeoTIFF data (clamped 0..1) if loaded, else fallback to random noise
            grid[idx].fuel = fuelLoaded ? std::clamp(rawFuel[idx], 0.0f, 1.0f)
                : 0.8f + dist(rng) * 0.2f;
        }
    }
    return grid;
}

// ----------------------------------------------------------------------------
// Compute Terrain Shading
// ----------------------------------------------------------------------------
// Maps grid elevation values to a normalized grayscale range [0..255] for renderer shading.
// Higher elevation cells map closer to 255 (white).
std::vector<std::uint8_t> computeTerrainShade(const std::vector<Cell>& grid) {
    float minElev = std::numeric_limits<float>::max();
    float maxElev = std::numeric_limits<float>::lowest();

    // Find global minimum and maximum elevation across the grid
    for (const Cell& c : grid) {
        minElev = std::min(minElev, c.elevation);
        maxElev = std::max(maxElev, c.elevation);
    }

    // Compute range with small epsilon to prevent division-by-zero on completely flat terrain
    const float elevRange = std::max(maxElev - minElev, 1e-6f);
    std::vector<std::uint8_t> shade(grid.size());

    // Normalize elevation to 8-bit unsigned byte intensity (0–255)
    for (size_t i = 0; i < grid.size(); ++i) {
        const float t = (grid[i].elevation - minElev) / elevRange;
        shade[i] = static_cast<std::uint8_t>(std::clamp(t, 0.0f, 1.0f) * 255.0f + 0.5f);
    }
    return shade;
}