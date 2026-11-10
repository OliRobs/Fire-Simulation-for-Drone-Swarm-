#include <algorithm>
#include <cmath>
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

// ============================================================================
// Main Application Entry Point
//
// Integrates simulation components:
// 1. LANDFIRE GeoTIFF ingestion (elevation, slope, FBFM40 fuel models, roads)
// 2. Rothermel surface fire spread on the cell grid
// 3. Distance Field Computation (Dynamic No-Fly Zones)
// 4. Fire Perimeter Tracking & Containment Boundaries
// 5. Autonomous Multi-Drone Mission Control & Path Re-Planning
// 6. Real-time Visualization & User Event Processing (SFML)
//
// Usage:  wildfire [--data <folder>] [--centre <lon> <lat>]
//   --data    folder holding the LANDFIRE downloads (searched recursively), default "landfire"
//   --centre  WGS84 longitude/latitude of the simulation window centre, default the raster centre
// Without LANDFIRE data it falls back to elevation_map.tif (+ optional fuel_map.tif) in the working folder.
// ============================================================================

namespace {
    // ------------------------------------------------------------------------
    // Nearest burnable cell to idx (spiral search), so the start fire doesn't land on a road or a town
    // ------------------------------------------------------------------------
    int nearestBurnable(const FireSimulation& sim, int idx) {
        const int cx = idx % GRID_WIDTH, cy = idx / GRID_WIDTH;
        for (int r = 0; r < GRID_WIDTH; ++r) {
            for (int y = cy - r; y <= cy + r; ++y) {
                for (int x = cx - r; x <= cx + r; ++x) {
                    if (std::max(std::abs(x - cx), std::abs(y - cy)) != r) continue;   // Ring of radius r only
                    if (x < 0 || x >= GRID_WIDTH || y < 0 || y >= GRID_HEIGHT) continue;
                    if (sim.canBurn(y * GRID_WIDTH + x)) return y * GRID_WIDTH + x;
                }
            }
        }
        return -1;   // Nothing on the map can burn (e.g. fuel far too wet)
    }

    // ------------------------------------------------------------------------
    // Print the weather assumptions: they are not in the LANDFIRE data, so make them visible
    // ------------------------------------------------------------------------
    void printWeather(const FireWeather& w) {
        std::printf("Fire weather (fixed scenario, not from data):\n"
                    "  20-ft wind %.0f km/h from %.0f deg, WAF %s\n"
                    "  fuel moisture: 1-h %.0f%%, 10-h %.0f%%, 100-h %.0f%%, live herb %.0f%%, live woody %.0f%%\n",
                    w.wind20ftKmh, w.windFromDeg,
                    w.windAdjustmentFactor > 0 ? std::to_string(w.windAdjustmentFactor).c_str() : "unsheltered (from fuel depth)",
                    w.moisture.dead1h * 100, w.moisture.dead10h * 100, w.moisture.dead100h * 100,
                    w.moisture.liveHerb * 100, w.moisture.liveWoody * 100);
    }
}

int main(int argc, char** argv) {
    // Initialize GDAL raster drivers for GeoTIFF file loading
    registerRasterDrivers();

    // ------------------------------------------------------------------------
    // Step 1: Load the LANDFIRE landscape and build the fire simulation
    // ------------------------------------------------------------------------
    std::string dataDir = "landfire";
    WindowRequest window;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--data") == 0 && i + 1 < argc) {
            dataDir = argv[++i];
        }
        else if ((std::strcmp(argv[i], "--centre") == 0 || std::strcmp(argv[i], "--center") == 0) && i + 2 < argc) {
            window.useLonLat = true;
            window.lon = std::atof(argv[++i]);
            window.lat = std::atof(argv[++i]);
        }
        else {
            std::cerr << "Usage: wildfire [--data <folder>] [--centre <lon> <lat>]\n";
            return 1;
        }
    }

    LandscapePaths paths = findLandfireLayers(dataDir);
    if (paths.elevation.empty() && paths.fuel.empty()) {
        // No LANDFIRE download: the original bundled DEM (and optional fuel map) in the working folder
        std::cout << "Falling back to elevation_map.tif / fuel_map.tif\n";
        paths.elevation = "elevation_map.tif";
        paths.fuel = "fuel_map.tif";
    }
    Terrain terrain = loadTerrain(paths, window);

    // Set up rendering engine with the elevation shading (white = highest)
    Renderer renderer(computeTerrainShade(terrain.grid));

    // Rothermel spread needs weather as well as the landscape; see config.h and README for what to supply
    const FireWeather weather;
    printWeather(weather);
    FireSimulation sim(std::move(terrain.grid), weather, terrain.cellSizeX, terrain.cellSizeY);

    // Initial full-grid visual repaint
    renderer.paintAll(sim.grid());

    // Seed the start fire at the centre of the map, or the nearest cell that can burn
    const int start = nearestBurnable(sim, (GRID_HEIGHT / 2) * GRID_WIDTH + GRID_WIDTH / 2);
    if (start >= 0) {
        sim.ignite(start);
        std::printf("Ignition at cell (%d, %d): fuel model %d, head fire %.1f m/min\n",
                    start % GRID_WIDTH, start / GRID_WIDTH, sim.grid()[start].fuelModel, sim.headRos(start));
    }
    else {
        std::cerr << "No burnable cell on the map (check the fuel layer and fuel moisture)\n";
    }

    // ------------------------------------------------------------------------
    // Step 2: Initialize Distance Field, Perimeter Tracking, and Mission Control
    // ------------------------------------------------------------------------
    // Live distance transform (EDT) maintaining safety buffers around active fire front
    DistanceTransform liveEdt(GRID_WIDTH, GRID_HEIGHT);

    // Manages committed containment perimeter boundaries
    Perimeter perimeter(sim.fireMask());

    // Central controller managing autonomous drone fleet, task allocation, and routing
    MissionController mission(liveEdt, perimeter.distance(),
        Vec2(2.5f, GRID_HEIGHT - 2.5f)); // Base station location (bottom-left)

    bool ringDirty = true; // Flag indicating UI/Visual perimeter ring needs updating
    int lastTitleMinute = -1; // Simulated minute last shown in the window title

    // ------------------------------------------------------------------------
    // Step 3: Main Execution & Simulation Event Loop
    // ------------------------------------------------------------------------
    while (renderer.window().isOpen()) {
        sf::Event event;

        // --- 3a. User Event Handling ---
        while (renderer.window().pollEvent(event)) {
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
                if (cx >= 0 && cx < GRID_WIDTH && cy >= 0 && cy < GRID_HEIGHT) {
                    sim.ignite(cy * GRID_WIDTH + cx);
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
            liveEdt.compute(sim.fireMask());

            // Advance containment perimeter sections if encroached upon
            if (perimeter.update(sim.fireMask())) {
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
        }

        // --- 3e. Render Frame ---
        renderer.draw(mission.plan(), mission.drones());

        // Show the simulated time in the title bar, refreshed once per simulated minute
        const int simMinute = static_cast<int>(sim.minutes());
        if (simMinute != lastTitleMinute) {
            lastTitleMinute = simMinute;
            char title[128];
            std::snprintf(title, sizeof(title), "Wildfire (Rothermel, LANDFIRE)  t = %dh %02dm  |  P perimeter, V mission, R roads",
                          simMinute / 60, simMinute % 60);
            renderer.window().setTitle(title);
        }
    }

    return 0;
}