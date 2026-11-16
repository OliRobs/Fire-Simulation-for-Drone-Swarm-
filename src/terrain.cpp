#include "terrain.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <cpl_vsi.h>
#include "fuel_models.h"
#include "raster_io.h"
#include "topography.h"

// ============================================================================
// Landscape Loader (LANDFIRE GeoTIFFs)
//
// Builds the fire environment from four LANDFIRE layers, all on one simulation grid:
//   Elevation (m)          -> Cell::elevation, and the aspect (upslope direction) via topography.cpp
//   Slope (degrees)        -> Cell::slopeDeg  (falls back to slope derived from the elevation)
//   FBFM40 fuel models     -> Cell::fuelModel (Scott & Burgan model number, looked up in fuel_models.cpp)
//   Operational roads      -> Cell::roadClass (fuel breaks, see ROADS_ARE_FUEL_BREAKS)
// plus, when present, USGS 3DEP 1 m DEM tiles, which replace the LANDFIRE elevation and slope: each 10 m cell's
// elevation, slope and aspect are worked out from the 100 one-metre pixels inside it.
// Missing layers fall back to simple defaults, with a message saying so.
// ============================================================================

namespace {
    namespace fs = std::filesystem;

    // LANDFIRE marks NoData with -9999 (and 32767 / 255 in some products). GDAL turns declared NoData into
    // NaN; this catches files where the value isn't declared.
    bool isNoData(float v) { return std::isnan(v) || v <= -9000.0f || v >= 32000.0f; }

    // Lower-case copy of a string, for case-insensitive file-name matching
    std::string lower(std::string s) {
        std::transform(s.begin(), s.end(), s.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
        return s;
    }

    // Where the DEM tile mosaic is kept (in GDAL's in-memory file system, so nothing is written to disk)
    const char* const DEM_MOSAIC_PATH = "/vsimem/wf_dem_mosaic.vrt";

    // Loads one layer onto the grid; returns false (and leaves `out` empty) if the path is empty or unusable
    bool loadLayer(const std::string& path, const RasterGrid& grid, Resample method, const char* what,
                   std::vector<float>& out) {
        out.clear();
        if (path.empty()) {
            std::cout << "No " << what << " layer found\n";
            return false;
        }
        return loadRaster(path.c_str(), grid, method, out);
    }
}

// ----------------------------------------------------------------------------
// Find the LANDFIRE layers in a folder
// ----------------------------------------------------------------------------
// Matches on the layer name inside the file name (e.g. LF2020_Elev_CONUS.tif, LF2025_FBFM40_CONUS.tif).
LandscapePaths findLandfireLayers(const std::string& dir, bool recursive) {
    LandscapePaths p;
    std::error_code ec;
    if (!fs::is_directory(dir, ec)) return p;

    // Matches one file against the layer names. "_CH_" and "_CBH_" are distinct, so the canopy layers can't
    // be confused with each other.
    auto consider = [&p](const fs::path& path) {
        std::string ext = path.extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char ch) { return static_cast<char>(std::tolower(ch)); });
        if (ext != ".tif" && ext != ".tiff") return;

        const std::string name = path.filename().string();
        const std::string full = path.string();
        if (name.find("_Elev_") != std::string::npos) p.elevation = full;
        else if (name.find("_SlpD_") != std::string::npos) p.slope = full;
        else if (name.find("FBFM40") != std::string::npos) p.fuel = full;
        else if (name.find("_Roads_") != std::string::npos) p.roads = full;
        else if (name.find("_CC_") != std::string::npos) p.canopyCover = full;
        else if (name.find("_CH_") != std::string::npos) p.canopyHeight = full;
        else if (name.find("_CBH_") != std::string::npos) p.canopyBaseHeight = full;
        else if (name.find("_CBD_") != std::string::npos) p.canopyBulkDensity = full;
    };
    if (recursive) {
        for (const fs::directory_entry& e : fs::recursive_directory_iterator(dir, ec))
            if (e.is_regular_file()) consider(e.path());
    }
    else {
        for (const fs::directory_entry& e : fs::directory_iterator(dir, ec))
            if (e.is_regular_file()) consider(e.path());
    }
    p.demTiles = findDemTiles(dir, recursive);
    if (!p.any()) return p;

    auto show = [](const char* label, const std::string& path) {
        std::cout << "  " << label << (path.empty() ? "(missing)" : path) << "\n";
    };
    std::cout << "LANDFIRE layers in '" << dir << "':\n";
    show("elevation:          ", p.elevation);
    show("slope:              ", p.slope);
    show("fuel (FBFM40):      ", p.fuel);
    show("roads:              ", p.roads);
    show("canopy cover:       ", p.canopyCover);
    show("canopy height:      ", p.canopyHeight);
    show("canopy base height: ", p.canopyBaseHeight);
    show("canopy bulk density:", p.canopyBulkDensity);
    std::cout << "  1 m DEM tiles:       " << p.demTiles.size() << "\n";
    for (const std::string& t : p.demTiles) std::cout << "    " << t << "\n";
    return p;
}

// ----------------------------------------------------------------------------
// Find the USGS 1 m DEM tiles in a folder
// ----------------------------------------------------------------------------
// USGS names its 1 m tiles USGS_1M_<UTM zone>_x<easting/10 km>y<northing/10 km>_<project>.tif, so any file
// with "USGS_1M" in its name is taken as a tile. Sorted so the mosaic is the same on every run.
std::vector<std::string> findDemTiles(const std::string& dir, bool recursive) {
    std::vector<std::string> tiles;
    std::error_code ec;
    if (!fs::is_directory(dir, ec)) return tiles;
    auto consider = [&tiles](const fs::path& path) {
        const std::string ext = lower(path.extension().string());
        if ((ext == ".tif" || ext == ".tiff") && lower(path.filename().string()).find("usgs_1m") != std::string::npos)
            tiles.push_back(path.string());
    };
    if (recursive) {
        for (const fs::directory_entry& e : fs::recursive_directory_iterator(dir, ec))
            if (e.is_regular_file()) consider(e.path());
    }
    else {
        for (const fs::directory_entry& e : fs::directory_iterator(dir, ec))
            if (e.is_regular_file()) consider(e.path());
    }
    std::sort(tiles.begin(), tiles.end());
    return tiles;
}

// ----------------------------------------------------------------------------
// Load Terrain Data
// ----------------------------------------------------------------------------
Terrain loadTerrain(const LandscapePaths& paths, const WindowRequest& window) {
    constexpr int N = GRID_WIDTH * GRID_HEIGHT;
    Terrain terrain;
    terrain.grid.resize(N);

    std::cout << "Loading GIS datasets via GDAL...\n";

    // --- 1. Simulation grid ---
    // With 1 m DEM tiles: the tiles' mosaic sets the grid (its UTM zone, centred on the tiles). Otherwise the
    // LANDFIRE elevation layer does (or the fuel layer if there's no DEM).
    RasterGrid full, grid;
    const bool haveDemTiles = !paths.demTiles.empty() && buildMosaic(paths.demTiles, DEM_MOSAIC_PATH);
    const std::string gridSource = haveDemTiles ? std::string(DEM_MOSAIC_PATH)
                                 : !paths.elevation.empty() ? paths.elevation : paths.fuel;
    bool haveGrid = !gridSource.empty() && readRasterGrid(gridSource.c_str(), GRID_WIDTH, GRID_HEIGHT, full);
    if (haveGrid) {
        grid = full;
        // Centre: requested lon/lat, else the middle of the raster
        double cx = 0.5 * (full.minX + full.maxX), cy = 0.5 * (full.minY + full.maxY);
        if (window.useLonLat && !lonLatToGrid(full, window.lon, window.lat, cx, cy)) {
            std::cerr << "Could not convert " << window.lon << ", " << window.lat
                      << " to the raster's coordinate system; using the raster centre\n";
            cx = 0.5 * (full.minX + full.maxX);
            cy = 0.5 * (full.minY + full.maxY);
        }
        int zone = 0;
        double rotation = 0.0;
        if (haveDemTiles && full.projected) {
            // The 1 m tiles are already north-up in UTM, so the window is cut straight out of them, its edges
            // on whole metres (four tiles x41-42 y379-380: easting 410-430 km, northing 3780-3800 km).
            // If the tiles are smaller than the window, it is still centred on them and the cells beyond the
            // DEM fall back to the LANDFIRE elevation and slope.
            if (!makeWindow(full, cx, cy, SIM_CELL_SIZE_M, grid)) {
                const double w = GRID_WIDTH * SIM_CELL_SIZE_M, h = GRID_HEIGHT * SIM_CELL_SIZE_M;
                grid.minX = std::round((cx - w / 2) / SIM_CELL_SIZE_M) * SIM_CELL_SIZE_M;
                grid.maxY = std::round((cy + h / 2) / SIM_CELL_SIZE_M) * SIM_CELL_SIZE_M;
                grid.maxX = grid.minX + w;
                grid.minY = grid.maxY - h;
                std::cout << "Warning: the DEM tiles are smaller than the " << w / 1000 << " km simulation window; "
                             "cells outside them use the LANDFIRE elevation\n";
            }
            std::printf("Simulation window: %d x %d cells of %.0f m on the 1 m DEM's grid, x %.0f .. %.0f, "
                        "y %.0f .. %.0f\n", GRID_WIDTH, GRID_HEIGHT, SIM_CELL_SIZE_M, grid.minX, grid.maxX,
                        grid.minY, grid.maxY);
        }
        // A large projected raster (LANDFIRE) gets a window of SIM_CELL_SIZE_M cells, north-up in the local UTM
        // zone (or on the file's own pixel grid if ALIGN_NORTH_UP is off); a small DEM is used whole
        else if (full.projected && ALIGN_NORTH_UP && makeNorthUpWindow(full, cx, cy, SIM_CELL_SIZE_M, grid, zone, rotation)) {
            std::printf("Simulation window: %d x %d cells of %.0f m, north-up on UTM zone %dN (the source grid is "
                        "rotated %.1f deg from north here, so layers are re-sampled onto the new grid)\n",
                        GRID_WIDTH, GRID_HEIGHT, SIM_CELL_SIZE_M, zone, rotation);
        }
        else if (full.projected && makeWindow(full, cx, cy, SIM_CELL_SIZE_M, grid)) {
            std::cout << "Simulation window: " << GRID_WIDTH << " x " << GRID_HEIGHT << " cells of "
                      << SIM_CELL_SIZE_M << " m on the source grid, x " << grid.minX << " .. " << grid.maxX
                      << ", y " << grid.minY << " .. " << grid.maxY << "\n";
        }
        else if (window.useLonLat) {
            std::cout << "Raster is smaller than the simulation window, so --centre is ignored\n";
        }
        terrain.cellSizeX = static_cast<float>(grid.cellSizeX());
        terrain.cellSizeY = static_cast<float>(grid.cellSizeY());
        std::cout << "Grid cell size: " << terrain.cellSizeX << " x " << terrain.cellSizeY << " (map units)\n";
        // The automaton draws square cells; a non-square footprint still works but looks stretched
        if (std::fabs(terrain.cellSizeX / terrain.cellSizeY - 1.0f) > 0.01f) {
            std::cerr << "Warning: grid cells are not square (stretched on screen)\n";
        }
    }

    // --- 2. Terrain from the 1 m DEM (if there are tiles) ---
    // Read now, so the LANDFIRE elevation and slope below are only needed for cells the tiles don't cover.
    std::vector<CellTopography> fineTopo;   // One entry per cell; empty without tiles
    double finePixel = 0.0;
    if (haveGrid && haveDemTiles && !loadFineTopography(DEM_MOSAIC_PATH, grid, fineTopo, finePixel)) {
        fineTopo.clear();
        std::cerr << "Could not use the 1 m DEM tiles; using the LANDFIRE elevation and slope instead\n";
    }
    auto fromFineDem = [&fineTopo](int i) { return !fineTopo.empty() && fineTopo[i].validPixels > 0; };
    int fineCells = 0;
    for (int i = 0; i < N; ++i) fineCells += fromFineDem(i) ? 1 : 0;

    // --- 3. Elevation ---
    // 1 m DEM mean where the tiles cover the cell; LANDFIRE elevation elsewhere (skipped if the tiles cover
    // the whole grid, which saves warping a 30 m layer that wouldn't be used).
    std::vector<float> raw;
    const bool needLandfireTerrain = fineCells < N;
    const bool haveLandfireElev = haveGrid && needLandfireTerrain
        && loadLayer(paths.elevation, grid, Resample::Average, "elevation", raw);
    if (haveLandfireElev) {
        // Cells with no data get the lowest valid elevation, so they don't create false cliffs
        float minValid = std::numeric_limits<float>::max();
        for (float v : raw) if (!isNoData(v)) minValid = std::min(minValid, v);
        for (int i = 0; i < N; ++i) terrain.grid[i].elevation = isNoData(raw[i]) ? minValid : raw[i];
    }
    else if (fineCells == 0) {
        std::cout << "No DEM: using a synthetic diagonal slope\n";
        for (int y = 0; y < GRID_HEIGHT; ++y)
            for (int x = 0; x < GRID_WIDTH; ++x)
                terrain.grid[y * GRID_WIDTH + x].elevation = (x * 0.3f) + ((GRID_HEIGHT - y) * 0.5f);
    }
    if (fineCells > 0) {
        // Cells the tiles miss (only if the tiles are smaller than the window) get the lowest DEM elevation
        // when there is no LANDFIRE elevation either
        float minFine = std::numeric_limits<float>::max();
        for (int i = 0; i < N; ++i) if (fromFineDem(i)) minFine = std::min(minFine, fineTopo[i].elevation);
        for (int i = 0; i < N; ++i) {
            if (fromFineDem(i)) terrain.grid[i].elevation = fineTopo[i].elevation;
            else if (!haveLandfireElev) terrain.grid[i].elevation = minFine;
        }
    }

    // --- 4. Slope and aspect ---
    // From the grid elevations first (Horn's method on the cells). Then, cell by cell, the best source wins:
    // the 1 m DEM (mean 1 m slope and overall uphill direction) where it covers the cell, else the measured
    // LANDFIRE slope (direction still from the elevations).
    computeTopography(terrain.grid, terrain.cellSizeX, terrain.cellSizeY, false);
    if (haveGrid && needLandfireTerrain && loadLayer(paths.slope, grid, Resample::Average, "slope", raw)) {
        for (int i = 0; i < N; ++i) {
            if (!fromFineDem(i) && !isNoData(raw[i]) && raw[i] >= 0.0f && raw[i] <= 90.0f) terrain.grid[i].slopeDeg = raw[i];
        }
    }
    else if (fineCells == 0) {
        std::cout << "Using slope derived from the elevation\n";
    }
    if (fineCells > 0) {
        double slopeSum = 0.0;
        for (int i = 0; i < N; ++i) {
            if (!fromFineDem(i)) continue;
            terrain.grid[i].slopeDeg = fineTopo[i].slopeDeg;
            terrain.grid[i].upslopeX = fineTopo[i].upslopeX;
            terrain.grid[i].upslopeY = fineTopo[i].upslopeY;
            slopeSum += fineTopo[i].slopeDeg;
        }
        std::printf("Terrain from the %.0f m DEM in %d of %d cells (mean slope %.1f deg)\n", finePixel, fineCells, N,
                    slopeSum / fineCells);
        fineTopo.clear();
        fineTopo.shrink_to_fit();   // ~80 MB for 4 million cells, no longer needed
    }
    VSIUnlink(DEM_MOSAIC_PATH);     // Free the in-memory mosaic description

    // --- 5. Fuel models (FBFM40 codes) ---
    if (haveGrid && loadLayer(paths.fuel, grid, Resample::Mode, "fuel", raw)) {
        std::map<int, int> counts;   // Fuel model -> number of cells, for the summary below
        for (int i = 0; i < N; ++i) {
            // Outside the fuel map: no fuel, so fire stops at the edge of the data
            const int code = isNoData(raw[i]) ? 0 : static_cast<int>(std::lround(raw[i]));
            terrain.grid[i].fuelModel = static_cast<std::int16_t>(code);
            ++counts[code];
        }
        // Report the most common fuel models so it's obvious what landscape was loaded
        std::vector<std::pair<int, int>> top(counts.begin(), counts.end());
        std::sort(top.begin(), top.end(), [](auto& a, auto& b) { return a.second > b.second; });
        std::cout << "Most common fuel models:";
        for (size_t k = 0; k < top.size() && k < 6; ++k) {
            const FuelModel* m = isNonburnableCode(top[k].first) ? nullptr : findFuelModel(top[k].first);
            std::cout << " " << (m ? m->name : (top[k].first <= 0 ? "none" : "NB")) << "(" << top[k].first << ") "
                      << (100 * top[k].second / N) << "%";
        }
        std::cout << "\n";
    }
    else {
        std::cout << "No fuel map: every cell gets fuel model " << DEFAULT_FUEL_MODEL << "\n";
        for (Cell& c : terrain.grid) c.fuelModel = DEFAULT_FUEL_MODEL;
    }

    // --- 6. Operational roads ---
    if (haveGrid && loadLayer(paths.roads, grid, Resample::Max, "roads", raw)) {
        int roadCells = 0;
        for (int i = 0; i < N; ++i) {
            const int v = isNoData(raw[i]) ? 0 : static_cast<int>(std::lround(raw[i]));
            const bool road = v >= ROAD_CLASS_MIN && v <= ROAD_CLASS_MAX;
            terrain.grid[i].roadClass = static_cast<std::uint8_t>(road ? v : 0);
            roadCells += road ? 1 : 0;
        }
        std::cout << "Road cells: " << roadCells << " (" << (ROADS_ARE_FUEL_BREAKS ? "fuel breaks" : "drawn only")
                  << ")\n";
    }

    // --- 7. Forest canopy (wind sheltering and crown fire) ---
    // LANDFIRE stores these as scaled integers and 0 for non-forest. Each layer is optional; a missing one
    // leaves the value at 0 ("no trees"), which switches off sheltering and crown fire for that cell.
    struct CanopyLayer {
        const std::string* path;
        const char* what;
        float scale;                 // Stored value x scale = real value
        float Cell::* field;
    };
    const CanopyLayer canopy[] = {
        { &paths.canopyCover, "canopy cover", 1.0f, &Cell::canopyCoverPct },          // percent
        { &paths.canopyHeight, "canopy height", 0.1f, &Cell::canopyHeightM },         // m x 10
        { &paths.canopyBaseHeight, "canopy base height", 0.1f, &Cell::canopyBaseHeightM }, // m x 10
        { &paths.canopyBulkDensity, "canopy bulk density", 0.01f, &Cell::canopyBulkDensity }, // kg/m^3 x 100
    };
    int canopyLayers = 0;
    for (const CanopyLayer& layer : canopy) {
        if (!(haveGrid && loadLayer(*layer.path, grid, Resample::Average, layer.what, raw))) continue;
        ++canopyLayers;
        for (int i = 0; i < N; ++i) {
            const float v = isNoData(raw[i]) ? 0.0f : std::max(raw[i], 0.0f);
            terrain.grid[i].*layer.field = v * layer.scale;
        }
    }
    if (canopyLayers > 0) {
        int forest = 0;
        for (const Cell& c : terrain.grid) forest += c.canopyCoverPct > 0.0f ? 1 : 0;
        std::cout << "Forest canopy over " << (100 * forest / N) << "% of the window\n";
        if (canopyLayers < 4) {
            std::cout << "Warning: only " << canopyLayers << " of the 4 canopy layers (CC, CH, CBH, CBD) found; "
                         "crown fire needs all of them\n";
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
