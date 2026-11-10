#pragma once
#include <string>
#include <vector>

// The only module that touches GDAL.
void registerRasterDrivers();

// Where the simulation grid sits on the ground: a coordinate system plus an extent, split into
// width x height cells. Every input layer is warped onto this one grid so the layers line up cell for cell.
struct RasterGrid {
    std::string wkt;                 // Coordinate system (WKT); empty if the file had none
    bool projected = false;          // True if the coordinate system is projected (units are metres, ft...)
    double minX = 0, minY = 0, maxX = 0, maxY = 0; // Extent in that coordinate system
    int width = 0, height = 0;       // Simulation grid size in cells
    // Pixel grid of the source file (used to line a window up with its pixels)
    double originX = 0, originY = 0; // Top-left corner of the file
    double pixelW = 1, pixelH = 1;   // Pixel size (both positive)

    double cellSizeX() const { return (maxX - minX) / width; }   // Ground size of one cell (metres if projected)
    double cellSizeY() const { return (maxY - minY) / height; }
};

// Reads a raster's coordinate system and full extent, to use as the simulation grid (normally the DEM's).
bool readRasterGrid(const char* filePath, int width, int height, RasterGrid& out);

// Cuts a width x height window of cellSize-sized cells out of `full`, centred as close to (centreX, centreY)
// as the file's extent allows, with its edges snapped to the file's pixel grid (so at native resolution
// every cell is exactly one source pixel). False if the file is smaller than the window.
bool makeWindow(const RasterGrid& full, double centreX, double centreY, double cellSize, RasterGrid& out);

// Like makeWindow, but the window is laid out on a NORTH-UP grid: the local UTM zone (NAD83, the LANDFIRE
// datum) chosen from the window centre. Needed for LANDFIRE because its CONUS Albers projection (EPSG:5070) is
// only north-up at 96 W; far from there its grid is rotated (13 degrees in Southern California), so the map,
// roads and wind directions came out tilted. Every layer is then warped onto this grid.
// (centreX, centreY) is given in full's coordinate system. The window is moved in from the edges of the file
// if needed. False if the file is too small or the coordinate systems can't be converted.
// rotationDeg returns how far the source grid was rotated from north at the centre (for the log).
bool makeNorthUpWindow(const RasterGrid& full, double centreX, double centreY, double cellSize, RasterGrid& out,
                       int& utmZone, double& rotationDeg);

// Converts WGS84 longitude/latitude (degrees) into the grid's coordinate system
bool lonLatToGrid(const RasterGrid& grid, double lon, double lat, double& x, double& y);

enum class Resample {
    Average,  // Continuous data (elevation): mean of the source pixels in each cell
    Mode,     // Categorical data (fuel-model codes): most common source value in each cell, never a blend
    Max       // Thin categorical features (roads): a cell takes the feature if any source pixel in it has one,
              // so a one-pixel road stays unbroken after the grid is rotated
};

// Reads band 1 of filePath, reprojected and resampled onto grid (one value per cell, row-major).
// Cells with no source data (outside the file, or NoData) come back as NaN.
// Returns false if the file can't be read or doesn't overlap the grid at all.
bool loadRaster(const char* filePath, const RasterGrid& grid, Resample method, std::vector<float>& outData);
