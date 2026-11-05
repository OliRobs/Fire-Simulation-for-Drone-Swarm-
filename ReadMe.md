# Wildfire simulator

Build (needs SFML 2.5+, GDAL, CMake 3.16+):

    cmake -S . -B build && cmake --build build
    ./build/wildfire        # run from the folder containing elevation_map.tif and fuel_map.tif
    ./build/wf_tests        # headless tests (no SFML / GDAL needed at runtime)

Keys: left-click ignites a cell, P toggles the perimeter ring, V toggles the drone mission overlay.

## Layout (src/)

| Module | Responsibility |
|---|---|
| config.h, types.h | constants, neighbour tables, Vec2, Cell |
| raster_io, terrain | GDAL GeoTIFF loading and grid/terrain-shade construction (only code using GDAL) |
| fire_model | wind/slope physics, FireSimulation (frontier list, precomputed ignition table) |
| distance_transform | exact Euclidean distance transform |
| perimeter | committed fire + sectional 20-cell perimeter updates |
| safe_planner | fire-safe legs: line of sight, A*, string pulling |
| voronoi, boustrophedon, aco | the three planning algorithms |
| mission_planner | zones, Voronoi regions, sweep sectors, ACO tour |
| drone, mission_controller | drone movement and safety; re-planning and route assignment |
| renderer | all SFML drawing |
| main.cpp | window, events, frame loop |

Library layering: wf_core (everything except GDAL/SFML) <- wf_io (GDAL), wf_render (SFML) <- wildfire.