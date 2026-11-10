# Wildfire simulator

Cellular-automaton wildfire on real terrain (GDAL GeoTIFFs), with a drone swarm that surveys round the fire:
boustrophedon sweeps inside the perimeter and an ant-colony-optimised tour outside it, both kept clear of
the fire.

## Build

Needs CMake 3.16+, **SFML 2.5 – 2.6** (not SFML 3, the API changed) and GDAL 3.x.

    cmake -S . -B build && cmake --build build
    ./build/wildfire        # run from the repo root (where elevation_map.tif is)
    ./build/wf_tests        # headless tests

To build just the simulation core and the tests, without SFML or GDAL:

    cmake -S . -B build -DWF_BUILD_APP=OFF && cmake --build build && ctest --test-dir build

Visual Studio: `UAVFireMap.vcxproj` expects GDAL in `C:\gdal` and SFML 2.6.1 in `C:\SFML-2.6.1`.

Keys: left-click ignites a cell, P toggles the perimeter ring, V toggles the drone mission overlay.

## Input data

| File | Used for | Notes |
|---|---|---|
| `elevation_map.tif` | Terrain and the simulation grid | Its coordinate system and extent define the grid (GRID_WIDTH x GRID_HEIGHT cells). The bundled file is British National Grid, 5 km square, so cells are 10 m. |
| `fuel_map.tif` (optional) | Fuel | Warped onto the DEM's grid, so a different projection or extent is fine as long as it overlaps. **FBFM40 codes** (LANDFIRE) are converted with the Scott & Burgan table in `fuel_models.cpp`; a map of values in 0..1 is used directly. Without it, fuel is random. |

LANDFIRE only covers the US, so it can't provide fuel for the bundled UK DEM. A fuel map for this area must come
from a UK source, converted to FBFM40 codes or 0..1 values.

## Layout (src/)

| Module | Responsibility |
|---|---|
| config.h, types.h | constants, neighbour tables, Vec2, Cell |
| raster_io, terrain | GDAL loading: warps every layer onto the DEM's grid (only code using GDAL) |
| fuel_models | FBFM40 code -> flammability and burn duration (Scott & Burgan 2005) |
| fire_model | wind/slope physics, FireSimulation (frontier list, precomputed ignition table) |
| distance_transform | exact Euclidean distance transform |
| perimeter | committed fire + sectional perimeter updates |
| safe_planner | fire-safe legs: line of sight, A*, string pulling; reports "no safe path" instead of guessing |
| voronoi, boustrophedon, aco | the three planning algorithms |
| mission_planner | zones, Voronoi regions, sweep sectors, ACO tour |
| drone, mission_controller | drone movement and safety; re-planning, re-routing blocked drones, route assignment |
| renderer | all SFML drawing |
| main.cpp | window, events, frame loop |

Library layering: wf_core (everything except GDAL/SFML) <- wf_io (GDAL), wf_render (SFML) <- wildfire.

## Key settings (config.h)

| Setting | Value | Meaning |
|---|---|---|
| PERIMETER_DISTANCE | 40 cells | Perimeter ring distance from the committed fire |
| PERIMETER_UPDATE_DISTANCE | 20 cells | A perimeter section moves once the fire gets this close to it |
| MIN_FIRE_DISTANCE | 10 cells | Survey standoff for sweeps and points of interest |
| HARD_MIN_FIRE_DISTANCE | 3 cells | Flight legs may cut inside the standoff, but never closer than this |
| SEARCH_DISTANCE | 70 cells | Outer edge of the ACO search area |

## Drone safety

- Every planned leg is checked against the fire; if a waypoint can only be reached by crossing the fire it is
  dropped from the route rather than flown in a straight line.
- Between re-plans the fire keeps growing, so each segment is re-checked as the drone reaches it. A blocked
  drone is re-routed straight away, or hovers if nothing is reachable, retrying every 15 frames.
- `keepClearOfFire` is the last-resort push away from the fire.
