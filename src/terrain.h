#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include "types.h"

// The simulation grid plus the ground size of one cell, which the slope and spread calculations need.
struct Terrain {
    std::vector<Cell> grid;
    float cellSizeX = DEFAULT_CELL_SIZE_M;   // Metres per cell, east-west
    float cellSizeY = DEFAULT_CELL_SIZE_M;   // Metres per cell, north-south
};

// GeoTIFF paths for the LANDFIRE layers. Any may be empty (missing): elevation falls back to a synthetic
// slope, slope is derived from elevation, fuel defaults to DEFAULT_FUEL_MODEL, roads to none and the canopy
// layers to "no trees" (no wind sheltering, no crown fire).
struct LandscapePaths {
    std::string elevation;  // LF20xx_Elev  : metres
    std::string slope;      // LF20xx_SlpD  : slope in degrees
    std::string fuel;       // LF20xx_FBFM40: Scott & Burgan fuel model numbers
    std::string roads;      // LF20xx_Roads : operational road classes 20-23
    std::string canopyCover;      // LF20xx_CC  : percent
    std::string canopyHeight;     // LF20xx_CH  : metres x 10
    std::string canopyBaseHeight; // LF20xx_CBH : metres x 10
    std::string canopyBulkDensity;// LF20xx_CBD : kg/m^3 x 100

    bool any() const { return !elevation.empty() || !fuel.empty(); }
};

// Searches `dir` for the layers by name (LFxxxx_<Layer>_CONUS.tif), so any LANDFIRE version year works.
// recursive: also look in sub-folders (the unzipped LANDFIRE download folders). Prints what it found.
LandscapePaths findLandfireLayers(const std::string& dir, bool recursive = true);

// Where to put the simulation window on a large raster
struct WindowRequest {
    bool useLonLat = false;   // false: centre of the elevation raster
    double lon = 0.0, lat = 0.0;  // WGS84 degrees
};

// Builds the grid from the GeoTIFFs. The elevation layer defines the coordinate system.
// Canopy layers are converted from LANDFIRE's stored integers to real units (RMRS-GTR-220: CBH is metres x 10,
// CBD is kg/m^3 x 100; LANDFIRE product pages: CH is metres x 10). If it is larger than
// GRID_WIDTH x GRID_HEIGHT cells of SIM_CELL_SIZE_M (LANDFIRE is), a window of that size is cut out at native
// resolution; a smaller DEM (such as the bundled elevation_map.tif) is used whole, as before. Every other layer
// is warped onto the same grid: average for elevation and slope, most common value for the category layers.
Terrain loadTerrain(const LandscapePaths& paths, const WindowRequest& window);

// Elevation mapped linearly to 0 (lowest, black) .. 255 (highest, white), one value per cell.
std::vector<std::uint8_t> computeTerrainShade(const std::vector<Cell>& grid);
