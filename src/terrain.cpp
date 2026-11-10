#include "terrain.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include "fuel_models.h"
#include "raster_io.h"

// ============================================================================
// Terrain and Elevation/Fuel Data Manager
//
// Responsible for loading GIS raster datasets (GeoTIFFs) into the 2D grid,
// applying fallback synthetic terrain generation if files are missing, and
// calculating the elevation shading used as the map background.
// ============================================================================

// ----------------------------------------------------------------------------
// Load Terrain Data
// ----------------------------------------------------------------------------
// 1. The DEM defines the simulation grid (its coordinate system and extent, split into GRID_WIDTH x GRID_HEIGHT).
// 2. Elevation is area-averaged onto that grid; the fuel map is warped onto the SAME grid using the most
//    common code per cell, so the two layers line up even if the fuel map has another projection or extent.
// 3. Anything missing falls back to procedural data, with a message saying so.
Terrain loadTerrain(const char* demPath, const char* fuelPath, std::mt19937& rng) {
    constexpr int N = GRID_WIDTH * GRID_HEIGHT;
    Terrain terrain;
    terrain.grid.resize(N);
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);

    std::cout << "Loading GIS datasets via GDAL...\n";

    // --- Simulation grid: the DEM's footprint, or the fuel map's if there is no DEM ---
    RasterGrid grid;
    bool haveGrid = readRasterGrid(demPath, GRID_WIDTH, GRID_HEIGHT, grid);
    if (!haveGrid) haveGrid = readRasterGrid(fuelPath, GRID_WIDTH, GRID_HEIGHT, grid);
    if (haveGrid) {
        terrain.cellSizeX = static_cast<float>(grid.cellSizeX());
        terrain.cellSizeY = static_cast<float>(grid.cellSizeY());
        std::cout << "Grid cell size: " << terrain.cellSizeX << " x " << terrain.cellSizeY << " (map units)\n";
    }

    // --- Elevation ---
    std::vector<float> rawElevation;
    const bool demLoaded = haveGrid && loadRaster(demPath, grid, Resample::Average, rawElevation);
    if (demLoaded) {
        // Cells with no data (NaN) get the lowest valid elevation, so they don't create false cliffs
        float minValid = std::numeric_limits<float>::max();
        for (float v : rawElevation) if (!std::isnan(v)) minValid = std::min(minValid, v);
        for (int i = 0; i < N; ++i) {
            terrain.grid[i].elevation = std::isnan(rawElevation[i]) ? minValid : rawElevation[i];
        }
    }
    else {
        std::cout << "No DEM: using a synthetic diagonal slope\n";
        for (int y = 0; y < GRID_HEIGHT; ++y)
            for (int x = 0; x < GRID_WIDTH; ++x)
                terrain.grid[y * GRID_WIDTH + x].elevation = (x * 0.3f) + ((GRID_HEIGHT - y) * 0.5f);
    }

    // --- Fuel ---
    std::vector<float> rawFuel;
    const bool fuelLoaded = haveGrid && loadRaster(fuelPath, grid, Resample::Mode, rawFuel);
    if (fuelLoaded) {
        // A map of FBFM40 codes has values like 91..204; a continuous fuel map stays within 0..1
        float maxValid = 0.0f;
        for (float v : rawFuel) if (!std::isnan(v)) maxValid = std::max(maxValid, v);
        const bool isFbfm40 = maxValid > 1.0f;
        std::cout << "Fuel map read as " << (isFbfm40 ? "FBFM40 fuel-model codes" : "0..1 flammability") << "\n";

        for (int i = 0; i < N; ++i) {
            Cell& c = terrain.grid[i];
            const float v = rawFuel[i];
            if (std::isnan(v)) {            // Outside the fuel map: treat as nonburnable
                c.flammability = 0.0f;
                c.fuel = 0.0f;
            }
            else if (isFbfm40) {            // Category code -> Scott & Burgan lookup
                const FuelProperties f = fuelFromFbfm40(static_cast<int>(std::lround(v)));
                c.flammability = f.flammability;
                c.fuel = f.load;
            }
            else {                          // Already a 0..1 value: use it for both
                c.flammability = std::clamp(v, 0.0f, 1.0f);
                c.fuel = c.flammability;
            }
        }
    }
    else {
        std::cout << "No usable fuel map: using random fuel (0.8 - 1.0)\n";
        for (Cell& c : terrain.grid) {
            c.fuel = 0.8f + dist(rng) * 0.2f;
            c.flammability = c.fuel;
        }
    }
    return terrain;
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