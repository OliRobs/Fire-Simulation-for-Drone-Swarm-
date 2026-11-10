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

// GeoTIFF paths for the four LANDFIRE layers. Any may be empty (missing): elevation falls back to a synthetic
// slope, slope is derived from elevation, fuel defaults to DEFAULT_FUEL_MODEL, and roads to none.
struct LandscapePaths {
    std::string elevation;  // LF20xx_Elev  : metres
    std::string slope;      // LF20xx_SlpD  : slope in degrees
    std::string fuel;       // LF20xx_FBFM40: Scott & Burgan fuel model numbers
    std::string roads;      // LF20xx_Roads : operational road classes 20-23
};

// Searches `dir` (and its sub-folders, as LANDFIRE downloads unzip to LFxxxx_<Layer>_CONUS/...) for the four
// layers by name, so any LANDFIRE version year works. Prints what it found.
LandscapePaths findLandfireLayers(const std::string& dir);

// Where to put the simulation window on a large raster
struct WindowRequest {
    bool useLonLat = false;   // false: centre of the elevation raster
    double lon = 0.0, lat = 0.0;  // WGS84 degrees
};

// Builds the grid from the GeoTIFFs. The elevation layer defines the coordinate system. If it is larger than
// GRID_WIDTH x GRID_HEIGHT cells of SIM_CELL_SIZE_M (LANDFIRE is), a window of that size is cut out at native
// resolution; a smaller DEM (such as the bundled elevation_map.tif) is used whole, as before. Every other layer
// is warped onto the same grid: average for elevation and slope, most common value for the category layers.
Terrain loadTerrain(const LandscapePaths& paths, const WindowRequest& window);

// Elevation mapped linearly to 0 (lowest, black) .. 255 (highest, white), one value per cell.
std::vector<std::uint8_t> computeTerrainShade(const std::vector<Cell>& grid);
