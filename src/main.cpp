#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include "distance_transform.h"
#include "fire_model.h"
#include "mission_controller.h"
#include "perimeter.h"
#include "raster_io.h"
#include "renderer.h"
#include "terrain.h"
#include "weather.h"
#include "wind_solver.h"

// ============================================================================
// Main Application Entry Point
//
// Integrates simulation components:
// 1. GeoTIFF ingestion: USGS 1 m DEM tiles (terrain) and LANDFIRE (fuel models, canopy, roads)
// 2. Terrain-adjusted wind (mass-conserving solver) and Rothermel surface fire spread on the cell grid
// 3. Distance Field Computation (Dynamic No-Fly Zones)
// 4. Fire Perimeter Tracking & Containment Boundaries
// 5. Autonomous Multi-Drone Mission Control & Path Re-Planning
// 6. Real-time Visualization & User Event Processing (SFML)
//
// Usage:  wildfire [--data <folder>] [--centre <lon> <lat>] [--weather <file>] [--wind <km/h>] [--wind-from <deg>]
//   --data       folder holding the LANDFIRE .tif files and the USGS 1 m DEM tiles (searched recursively).
//                Default: ./landfire, then the working folder itself (files such as LF2020_Elev_CONUS.tif placed
//                in the repo root); DEM tiles (USGS_1M_*.tif) are also looked for in ./dem
//   --centre     WGS84 longitude/latitude of the simulation window centre, default the centre of the DEM tiles
//                (or of the LANDFIRE raster without them)
//   --weather    weather scenario file (key = value lines, see scenarios/default.ini and weather.h); keys it
//                leaves out keep config.h's defaults
//   --wind       open 20-ft wind speed, km/h, and --wind-from its direction (blowing FROM, degrees clockwise
//                from north); applied after --weather
// Without LANDFIRE data it falls back to elevation_map.tif (+ optional fuel_map.tif) in the working folder.
// ============================================================================

namespace {
    // ------------------------------------------------------------------------
    // Nearest burnable cell to idx (spiral search), so the start fire doesn't land on a road or a town
    // ------------------------------------------------------------------------
    int nearestBurnable(const FireSimulation& sim, int idx) {
        const GridSpec& g = sim.gridSpec();
        const int cx = g.xOf(idx), cy = g.yOf(idx);
        for (int r = 0; r < std::max(g.width, g.height); ++r) {
            for (int y = cy - r; y <= cy + r; ++y) {
                for (int x = cx - r; x <= cx + r; ++x) {
                    if (std::max(std::abs(x - cx), std::abs(y - cy)) != r) continue;   // Ring of radius r only
                    if (!g.inside(x, y)) continue;
                    if (sim.canBurn(g.index(x, y))) return g.index(x, y);
                }
            }
        }
        return -1;   // Nothing on the map can burn (e.g. fuel far too wet)
    }

    // ------------------------------------------------------------------------
    // Print the weather assumptions: they are not in the LANDFIRE data, so make them visible
    // ------------------------------------------------------------------------
    // source: where the values came from (the scenario file, or config.h's defaults), plus any overrides
    void printWeather(const FireWeather& w, const std::string& source) {
        std::printf("Fire weather (fixed scenario, not from data; %s):\n"
                    "  20-ft wind %.0f km/h from %.0f deg, WAF %s\n"
                    "  fuel moisture: 1-h %.0f%%, 10-h %.0f%%, 100-h %.0f%%, live herb %.0f%%, live woody %.0f%%\n"
                    "  foliar moisture %.0f%%, crown fire %s, wind limit %s, terrain wind stability alpha %.2f\n",
                    source.c_str(), w.wind20ftKmh, w.windFromDeg,
                    w.windAdjustmentFactor > 0 ? std::to_string(w.windAdjustmentFactor).c_str() : "from fuel depth and canopy (Andrews 2012)",
                    w.moisture.dead1h * 100, w.moisture.dead10h * 100, w.moisture.dead100h * 100,
                    w.moisture.liveHerb * 100, w.moisture.liveWoody * 100,
                    w.foliarMoisturePct, w.crownFire ? "on" : "off", w.windLimit ? "on" : "off",
                    w.windStabilityAlpha);
    }
}

int main(int argc, char** argv) {
    // Initialize GDAL raster drivers for GeoTIFF file loading
    registerRasterDrivers();

    // ------------------------------------------------------------------------
    // Step 1: Load the LANDFIRE landscape and build the fire simulation
    // ------------------------------------------------------------------------
    std::string dataDir;   // Empty: look in ./landfire, then next to the program's working folder
    WindowRequest window;
    std::string weatherFile;                  // Empty: config.h's default scenario
    double windOverride = -1.0, windFromOverride = -1.0;   // < 0: not given
    // A command-line number, or exit with a message if it isn't one (strtod must use the whole argument)
    auto numberArg = [](const char* option, const char* text) {
        char* end = nullptr;
        const double v = std::strtod(text, &end);
        if (end == text || *end != '\0' || !std::isfinite(v) || v < 0.0) {
            std::cerr << option << " needs a number >= 0, got '" << text << "'\n";
            std::exit(1);
        }
        return v;
    };
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--data") == 0 && i + 1 < argc) {
            dataDir = argv[++i];
        }
        else if ((std::strcmp(argv[i], "--centre") == 0 || std::strcmp(argv[i], "--center") == 0) && i + 2 < argc) {
            window.useLonLat = true;
            window.lon = std::atof(argv[++i]);
            window.lat = std::atof(argv[++i]);
        }
        else if (std::strcmp(argv[i], "--weather") == 0 && i + 1 < argc) {
            weatherFile = argv[++i];
        }
        else if (std::strcmp(argv[i], "--wind") == 0 && i + 1 < argc) {
            windOverride = numberArg("--wind", argv[++i]);
        }
        else if (std::strcmp(argv[i], "--wind-from") == 0 && i + 1 < argc) {
            windFromOverride = numberArg("--wind-from", argv[++i]);
            if (windFromOverride > 360.0) { std::cerr << "--wind-from needs 0 .. 360 degrees\n"; return 1; }
        }
        else {
            // --help / -h: the usage on stdout, success; anything else: the usage on stderr, failure
            const bool help = std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0;
            (help ? std::cout : std::cerr)
                << "Usage: wildfire [--data <folder>] [--centre <lon> <lat>] [--weather <file>] [--wind <km/h>]"
                   " [--wind-from <deg>]\n"
                   "Weather file keys (key = value, # comments; see scenarios/default.ini):\n" << weatherKeysHelp();
            return help ? 0 : 1;
        }
    }

    // --- Weather scenario: config.h's defaults, then the file, then the command-line overrides ---
    FireWeather weather = defaultWeather();
    std::string weatherSource = "config.h defaults";
    if (!weatherFile.empty()) {
        try {
            weather = loadWeatherFile(weatherFile);
        }
        catch (const WeatherFileError& e) {
            std::cerr << "Weather file error: " << e.what() << "\n";
            return 1;
        }
        weatherSource = "from " + weatherFile;
    }
    if (windOverride >= 0.0) {
        weather.wind20ftKmh = static_cast<float>(windOverride);
        weatherSource += ", --wind";
    }
    if (windFromOverride >= 0.0) {
        weather.windFromDeg = static_cast<float>(windFromOverride);
        weatherSource += ", --wind-from";
    }

    // Where the LANDFIRE .tif files are: --data folder, else a landfire/ folder (and its sub-folders), else
    // loose in the working folder (the repo root when run from Visual Studio or as ./build/wildfire)
    LandscapePaths paths;
    if (!dataDir.empty()) {
        paths = findLandfireLayers(dataDir, true);
        if (!paths.any()) std::cerr << "No LANDFIRE layers found in '" << dataDir << "'\n";
    }
    else {
        paths = findLandfireLayers("landfire", true);
        if (!paths.any()) paths = findLandfireLayers(".", false);
    }
    // The USGS 1 m DEM tiles may sit apart from the LANDFIRE layers: also look in a dem/ folder and loose in
    // the working folder (they are git-ignored there too)
    if (paths.demTiles.empty() && dataDir.empty()) {
        paths.demTiles = findDemTiles("dem", true);
        if (paths.demTiles.empty()) paths.demTiles = findDemTiles(".", false);
        if (!paths.demTiles.empty()) std::cout << "1 m DEM tiles found: " << paths.demTiles.size() << "\n";
    }
    if (!paths.any()) {
        // No LANDFIRE download: the original bundled DEM (and optional fuel map) in the working folder
        std::cout << "Falling back to elevation_map.tif / fuel_map.tif\n";
        paths.elevation = "elevation_map.tif";
        paths.fuel = "fuel_map.tif";
    }
    // The app's grid size is config.h's GRID_WIDTH x GRID_HEIGHT; everything below takes it from terrain.spec
    Terrain terrain = loadTerrain(paths, window, GRID_WIDTH, GRID_HEIGHT);
    const GridSpec grid = terrain.spec;

    // Set up rendering engine with the elevation shading (white = highest)
    Renderer renderer(grid, computeTerrainShade(terrain.grid));

    // Rothermel spread needs weather as well as the landscape (read above; see weather.h and README)
    printWeather(weather, weatherSource);

    // Terrain-adjusted wind: the open wind is adjusted over the terrain once, before the fire starts, so each
    // cell gets its own 20-ft wind speed and direction (wind_solver.h). Off: one wind for the whole map.
    WindField windField;
    if (USE_TERRAIN_WIND) {
        const auto t0 = std::chrono::steady_clock::now();
        WindSolverSettings windSettings;
        windSettings.alpha = weather.windStabilityAlpha;
        windField = solveTerrainWind(grid, terrain.grid, weather.wind20ftKmh, weather.windFromDeg, windSettings);
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        const auto range = std::minmax_element(windField.speedKmh.begin(), windField.speedKmh.end());
        std::printf("Terrain wind: %d x %d x %d mesh, %d CG iterations (residual %.1e%s), %.1f s\n"
                    "  20-ft wind %.1f to %.1f km/h across the map; mass error (RMS div) %.2e -> %.2e 1/s\n",
                    windField.columnsX, windField.columnsY, windField.layers, windField.iterations,
                    windField.relResidual, windField.converged ? "" : ", NOT converged", seconds,
                    *range.first, *range.second, windField.divergenceBefore, windField.divergenceAfter);
    }
    FireSimulation sim(grid, std::move(terrain.grid), weather, SIM_MINUTES_PER_TICK,
                       USE_TERRAIN_WIND ? &windField : nullptr);

    // The wind arrows ('W') draw this field; their colours compare each cell with the open wind
    renderer.setWindField(USE_TERRAIN_WIND ? &windField : nullptr, weather.wind20ftKmh);

    // Initial full-grid visual repaint
    renderer.paintAll(sim.grid());

    // Seed the start fire at the centre of the map, or the nearest cell that can burn
    const int start = nearestBurnable(sim, grid.index(grid.width / 2, grid.height / 2));
    // Summarise how the landscape would burn under this weather (crown fire needs the canopy layers)
    int passive = 0, active = 0, burnableCells = 0;
    for (int i = 0; i < static_cast<int>(grid.count()); ++i) {
        if (!sim.canBurn(i)) continue;
        ++burnableCells;
        passive += sim.grid()[i].fireType == static_cast<std::uint8_t>(FireType::Passive) ? 1 : 0;
        active += sim.grid()[i].fireType == static_cast<std::uint8_t>(FireType::Active) ? 1 : 0;
    }
    std::printf("Burnable cells: %d. Crown fire potential: %d torching (passive), %d active crown fire\n",
                burnableCells, passive, active);

    if (start >= 0) {
        sim.ignite(start);
        std::printf("Ignition at cell (%d, %d): fuel model %d, head fire %.1f m/min\n",
                    grid.xOf(start), grid.yOf(start), sim.grid()[start].fuelModel, sim.headRos(start));
    }
    else {
        std::cerr << "No burnable cell on the map (check the fuel layer and fuel moisture)\n";
    }

    // ------------------------------------------------------------------------
    // Step 2: Initialize Distance Field, Perimeter Tracking, and Mission Control
    // ------------------------------------------------------------------------
    // Live distance transform (EDT) maintaining safety buffers around active fire front
    // Only worked out within LIVE_EDT_RANGE of the fire: that is as far as any reader looks. Incremental: each
    // frame only the cells that caught fire since the last one are added (the fire never shrinks), and distances
    // beyond the range read "very far" (distance_transform.h). The readers (SafePlanner, keepClearOfFire, the
    // planner's zone test, the renderer's no-fly holes) compare against NO_FLY_SQ or HARD_NO_FLY_SQ, which the
    // static_assert below keeps inside the range.
    static_assert(NO_FLY_SQ < (LIVE_EDT_RANGE + 1.0) * (LIVE_EDT_RANGE + 1.0) && HARD_NO_FLY_SQ <= NO_FLY_SQ,
                  "the live distance field must reach every threshold its readers compare it with");
    DistanceTransform liveEdt(grid, LIVE_EDT_RANGE, DistanceTransform::Mode::IncrementalClamped);

    // Manages committed containment perimeter boundaries
    Perimeter perimeter(grid, sim.fireMask());

    // Central controller managing autonomous drone fleet, task allocation, and routing
    MissionController mission(liveEdt, perimeter.distance(),
        Vec2(2.5f, grid.height - 2.5f)); // Base station location (bottom-left)

    bool ringDirty = true; // Flag indicating UI/Visual perimeter ring needs updating
    int lastTitleMinute = -1; // Simulated minute last shown in the window title
    float lastTitleZoom = 0.0f; // Zoom last shown in the window title

    // ------------------------------------------------------------------------
    // Step 3: Main Execution & Simulation Event Loop
    // ------------------------------------------------------------------------
    while (renderer.window().isOpen()) {
        sf::Event event;

        // --- 3a. User Event Handling ---
        while (renderer.window().pollEvent(event)) {
            // Zoom / pan controls (mouse wheel, right-drag, +/-, arrows, 0) are handled by the renderer
            if (renderer.handleViewEvent(event)) continue;

            if (event.type == sf::Event::Closed) {
                renderer.window().close();
            }

            // Left-Click: Manually ignite a new spot fire at mouse position.
            // mapPixelToCoords converts window pixels to world coordinates, so this stays correct
            // after the window has been resized.
            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
                const sf::Vector2f world = renderer.window().mapPixelToCoords(
                    sf::Vector2i(event.mouseButton.x, event.mouseButton.y));
                const int cx = static_cast<int>(std::floor(world.x / CELL_SIZE));
                const int cy = static_cast<int>(std::floor(world.y / CELL_SIZE));
                if (grid.inside(cx, cy)) {
                    sim.ignite(grid.index(cx, cy));
                }
            }

            // 'P' Key: Toggle visual rendering of containment perimeter ring
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::P) {
                renderer.showPerimeter = !renderer.showPerimeter;
                ringDirty = true;
            }

            // 'V' Key: Toggle visual overlay of drone mission routes and waypoints
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::V) {
                renderer.showMission = !renderer.showMission;
            }

            // 'F' Key: Zoom to fit the fire (burning and burned cells)
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::F) {
                renderer.frameMask(sim.fireMask());
            }

            // 'W' Key: Toggle the terrain wind arrows (only when a terrain wind field was solved)
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::W) {
                if (!renderer.hasWindField()) {
                    std::printf("No terrain wind field to show (USE_TERRAIN_WIND is off in config.h)\n");
                }
                else {
                    renderer.showWind = !renderer.showWind;
                    if (renderer.showWind) {
                        std::printf("Wind arrows: 20-ft wind, length by speed (longest = %.1f km/h); white = the open "
                                    "%.0f km/h wind, blue = slower, magenta = faster (full colour at 30 %%)\n",
                                    renderer.windMaxKmh(), weather.wind20ftKmh);
                    }
                }
            }

            // 'R' Key: Toggle the road overlay (repaints the base map, then the perimeter ring on top)
            if (event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::R) {
                renderer.showRoads = !renderer.showRoads;
                renderer.paintAll(sim.grid());
                ringDirty = true;
            }
        }

        // --- 3b. Advance Wildfire Simulation ---
        sim.step();

        // Incremental visual render update for cells whose burning state changed
        renderer.paintCells(sim.grid(), sim.takeDirtyCells());

        // --- 3c. Handle Dynamic Fire Growth & Spatial Updates ---
        if (sim.takeFireGrew()) {
            // Recalculate dynamic Euclidean distance field for no-fly safety zones
            // (fireBounds() tells it where the fire is, so it needn't scan the grid to find it)
            liveEdt.compute(sim.fireMask(), sim.fireBounds());

            // Advance containment perimeter sections if encroached upon (only looks inside the fire's box)
            if (perimeter.update(sim.fireMask(), sim.fireBounds())) {
                ringDirty = true;
            }

            // Notify mission controller to trigger dynamic path re-planning
            mission.markFireChanged();
        }

        // Update containment perimeter visual ring overlay if modified
        if (ringDirty) {
            renderer.setRing(perimeter.ringCells(), sim.grid());
            ringDirty = false;
        }

        // --- 3d. Update Mission Controller & Drone Navigation ---
        mission.update();

        // Update mission graphics (routes, waypoints, no-fly bounds) if mission plan changed
        if (mission.takePlanChanged()) {
            renderer.setMissionVisuals(mission.plan(), mission.drones(), liveEdt);
            lastTitleMinute = -1;   // Refresh the title so the new drone counts show straight away
        }

        // --- 3e. Render Frame ---
        renderer.draw(mission.plan(), mission.drones());

        // Show the simulated time, zoom and drone counts in the title bar, refreshed each simulated minute, on
        // zooming, or after a re-plan
        const int simMinute = static_cast<int>(sim.minutes());
        if (simMinute != lastTitleMinute || renderer.zoom() != lastTitleZoom) {
            lastTitleMinute = simMinute;
            lastTitleZoom = renderer.zoom();
            char title[192];
            std::snprintf(title, sizeof(title),
                          "Wildfire  t = %dh %02dm  zoom %.1fx  |  wheel zoom, right-drag pan, F fire, 0 reset, P V R W overlays",
                          simMinute / 60, simMinute % 60, renderer.zoom());
            // Followed by the drone counts (and any over-limit warning) from the latest plan
            renderer.window().setTitle(std::string(title) + renderer.fleetStatus());
        }
    }

    return 0;
}