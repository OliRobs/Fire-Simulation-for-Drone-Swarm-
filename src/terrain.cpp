#include "terrain.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
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
// Missing layers fall back to simple defaults, with a message saying so.
// ============================================================================

namespace {
    namespace fs = std::filesystem;

    // LANDFIRE marks NoData with -9999 (and 32767 / 255 in some products). GDAL turns declared NoData into
    // NaN; this catches files where the value isn't declared.
    bool isNoData(float v) { return std::isnan(v) || v <= -9000.0f || v >= 32000.0f; }

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
    return p;
}

// ----------------------------------------------------------------------------
// Load Terrain Data
// ----------------------------------------------------------------------------
Terrain loadTerrain(const LandscapePaths& paths, const WindowRequest& window) {
    constexpr int N = GRID_WIDTH * GRID_HEIGHT;
    Terrain terrain;
    terrain.grid.resize(N);

    std::cout << "Loading GIS datasets via GDAL...\n";

    // --- 1. Simulation grid: from the elevation layer (or the fuel layer if there's no DEM) ---
    RasterGrid full, grid;
    const std::string& gridSource = !paths.elevation.empty() ? paths.elevation : paths.fuel;
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
        // A large projected raster (LANDFIRE) gets a native-resolution window; a small DEM is used whole
        if (full.projected && makeWindow(full, cx, cy, SIM_CELL_SIZE_M, grid)) {
            std::cout << "Simulation window: " << GRID_WIDTH << " x " << GRID_HEIGHT << " cells of "
                      << SIM_CELL_SIZE_M << " m, x " << grid.minX << " .. " << grid.maxX
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

    // --- 2. Elevation ---
    std::vector<float> raw;
    if (haveGrid && loadLayer(paths.elevation, grid, Resample::Average, "elevation", raw)) {
        // Cells with no data get the lowest valid elevation, so they don't create false cliffs
        float minValid = std::numeric_limits<float>::max();
        for (float v : raw) if (!isNoData(v)) minValid = std::min(minValid, v);
        for (int i = 0; i < N; ++i) terrain.grid[i].elevation = isNoData(raw[i]) ? minValid : raw[i];
    }
    else {
        std::cout << "No DEM: using a synthetic diagonal slope\n";
        for (int y = 0; y < GRID_HEIGHT; ++y)
            for (int x = 0; x < GRID_WIDTH; ++x)
                terrain.grid[y * GRID_WIDTH + x].elevation = (x * 0.3f) + ((GRID_HEIGHT - y) * 0.5f);
    }

    // --- 3. Slope and aspect ---
    // Aspect (which way is uphill) always comes from the DEM. Slope steepness is derived from it first, then
    // replaced by the measured LANDFIRE slope wherever that layer has data.
    computeTopography(terrain.grid, terrain.cellSizeX, terrain.cellSizeY, false);
    if (haveGrid && loadLayer(paths.slope, grid, Resample::Average, "slope", raw)) {
        for (int i = 0; i < N; ++i) {
            if (!isNoData(raw[i]) && raw[i] >= 0.0f && raw[i] <= 90.0f) terrain.grid[i].slopeDeg = raw[i];
        }
    }
    else {
        std::cout << "Using slope derived from the elevation\n";
    }

    // --- 4. Fuel models (FBFM40 codes) ---
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

    // --- 5. Operational roads ---
    if (haveGrid && loadLayer(paths.roads, grid, Resample::Mode, "roads", raw)) {
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

    // --- 6. Forest canopy (wind sheltering and crown fire) ---
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
