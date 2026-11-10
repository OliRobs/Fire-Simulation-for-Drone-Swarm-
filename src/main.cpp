#include <cmath>
#include <random>
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
// 1. Terrain & GeoTIFF Raster Data Ingestion
// 2. Stochastic Wildfire Simulation Loop
// 3. Distance Field Computation (Dynamic No-Fly Zones)
// 4. Fire Perimeter Tracking & Containment Boundaries
// 5. Autonomous Multi-Drone Mission Control & Path Re-Planning
// 6. Real-time Visualization & User Event Processing (SFML)
// ============================================================================
int main() {
    // Initialize GDAL raster drivers for GeoTIFF file loading
    registerRasterDrivers();

    // ------------------------------------------------------------------------
    // Step 1: Initialize Terrain, Shading, and Wildfire Simulation
    // ------------------------------------------------------------------------
    std::mt19937 rng(1337);

    // Load real-world elevation and fuel datasets into the cell grid (the DEM defines the grid)
    Terrain terrain = loadTerrain("elevation_map.tif", "fuel_map.tif", rng);

    // Set up rendering engine with the elevation shading (white = highest)
    Renderer renderer(computeTerrainShade(terrain.grid));

    // Instantiate wildfire engine with vector wind (blowing North-East: x=0.707, y=-0.707).
    // The cell size lets the slope term use real rise-over-run.
    FireSimulation sim(std::move(terrain.grid), Vec2(0.707f, -0.707f), 1337,
                       terrain.cellSizeX, terrain.cellSizeY);

    // Initial full-grid visual repaint
    renderer.paintAll(sim.grid());

    // Seed initial ignition point at the exact center of the grid map
    sim.ignite((GRID_HEIGHT / 2) * GRID_WIDTH + GRID_WIDTH / 2);

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
    }

    return 0;
}