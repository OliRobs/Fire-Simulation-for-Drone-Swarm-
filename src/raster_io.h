#pragma once
#include <string>
#include <vector>

// The only module that touches GDAL.
void registerRasterDrivers();

// Where the simulation grid sits on the ground: a coordinate system plus an extent, split into
// width x height cells. Every input layer is warped onto this one grid so the layers line up cell for cell.
struct RasterGrid {
    std::string wkt;                 // Coordinate system (WKT); empty if the file had none
    double minX = 0, minY = 0, maxX = 0, maxY = 0; // Extent in that coordinate system
    int width = 0, height = 0;       // Simulation grid size in cells

    double cellSizeX() const { return (maxX - minX) / width; }   // Ground size of one cell (metres if projected)
    double cellSizeY() const { return (maxY - minY) / height; }
};

// Reads a raster's coordinate system and full extent, to use as the simulation grid (normally the DEM's).
bool readRasterGrid(const char* filePath, int width, int height, RasterGrid& out);

enum class Resample {
    Average,  // Continuous data (elevation): mean of the source pixels in each cell
    Mode      // Categorical data (fuel-model codes): most common source value in each cell, never a blend
};

// Reads band 1 of filePath, reprojected and resampled onto grid (one value per cell, row-major).
// Cells with no source data (outside the file, or NoData) come back as NaN.
// Returns false if the file can't be read or doesn't overlap the grid at all.
bool loadRaster(const char* filePath, const RasterGrid& grid, Resample method, std::vector<float>& outData);
