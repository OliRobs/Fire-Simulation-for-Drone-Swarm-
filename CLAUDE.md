# Wildfire simulator with drone survey

C++17 cellular-automaton wildfire on real terrain (GDAL GeoTIFFs, drawn with SFML), plus a drone swarm that
surveys round the fire. README.md has the module table and input-data notes.

## Build and test

- Full build: `cmake -S . -B build && cmake --build build`, then run `./build/wildfire` from the repo root
  (it loads the LANDFIRE layers from the repo root or `landfire/`, falling back to `elevation_map.tif`).
- Tests: `./build/wf_tests` (or `ctest --test-dir build`). They are headless and take about 5 seconds.
- Core and tests only, without SFML or GDAL: add `-DWF_BUILD_APP=OFF`.
- The owner builds day to day with Visual Studio: `UAVFireMap.vcxproj`, SFML in `C:\SFML-2.6.1`, GDAL in
  `C:\gdal`. The project is x64 only (the SFML and GDAL builds are 64-bit). Debugger settings (working
  directory, PATH to the DLLs) live in the .vcxproj itself, because `.vcxproj.user` is not tracked.
  Adding or removing a source file means updating `CMakeLists.txt`, `UAVFireMap.vcxproj` AND
  `UAVFireMap.vcxproj.filters`.
- CI (`.github/workflows/build.yml`): full build and tests on Ubuntu 24.04, plus an MSVC build of the core and
  tests on Windows. Both must pass.

## Constraints

- SFML stays on 2.6. The code uses the SFML 2 API (`sf::Uint8`, `VideoMode(w, h)`, `pollEvent(event)`);
  don't port to SFML 3, which would break the Visual Studio setup. vcpkg's `sfml` port is SFML 3, so don't
  use it.
- Library layering: only `raster_io` and `terrain` use GDAL; only `renderer` and `main.cpp` use SFML.
  Everything else is in `wf_core`, which must build without either so the tests stay headless.
- The values in `config.h` (perimeter distances, fire standoffs, drone counts, speeds) are design decisions.
  Ask before changing them.
- Input data is LANDFIRE (US only, EPSG:5070 Albers, 30 m): elevation, slope (degrees), FBFM40 fuel models,
  operational roads and the canopy layers CC (%), CH (m x 10), CBH (m x 10), CBD (kg/m3 x 100). Files are found by
  name in `--data`, else `landfire/`, else loose in the working folder (the repo root). The elevation layer defines the simulation grid: a
  GRID_WIDTH x GRID_HEIGHT window at native 30 m, snapped to the pixel grid; every other layer is warped onto it.
  The .tif files are git-ignored. `elevation_map.tif` (UK, EPSG:27700, 5 km) is only the no-data fallback.
- Fuel rasters hold FBFM40 category codes: look them up through `fuel_models.cpp`, never average or clamp them
  (resample with Mode). Same for road classes.
- Fire spread is Rothermel (`rothermel.cpp`), equations from Andrews 2018 (RMRS-GTR-371); canopy wind sheltering
  from Andrews 2012; crown fire (`crown_fire.cpp`) is the Finney 1998 method (Van Wagner 1977 thresholds, Rothermel
  1991 crown rate). Weather (fuel moisture, foliar moisture, wind) is not in the LANDFIRE data; the placeholders
  are in `config.h` / `FireWeather` (README lists them).

## Drone survey design

These behaviours were asked for explicitly. Keep them when changing planning code:
- Drones must never fly over the fire. `SafePlanner::leg` returns an empty path when there is no safe route;
  callers skip the waypoint or hold position, never fall back to a straight line.
- Legs may cut inside the `MIN_FIRE_DISTANCE` survey standoff when the path allows, but never closer than
  `HARD_MIN_FIRE_DISTANCE`.
- Boustrophedon sweeps cover inside the perimeter; the ant-colony (ACO) tour covers from the perimeter out to
  `SEARCH_DISTANCE`. The search area is split into Voronoi regions that give the points of interest.
- ACO drones first fly to the closest point on the tour, then follow it, evenly spaced and in one direction.
- The perimeter updates in sections, only where the fire front has come within `PERIMETER_UPDATE_DISTANCE`.

## Conventions

- Comment all code clearly: a header block above each function saying what it does and why, and comments on
  non-obvious steps. Match the existing style (`// ----` section banners).
- Any change to fire spread, planning or drone movement needs a test in `tests/tests.cpp`. Run `wf_tests`
  before committing.
- Work on a branch and open a pull request; don't push to `main`.
- Physics and fuel constants come from the papers in the claude.ai project (Scott & Burgan 2005 for the fuel
  models; Andrews 2018 RMRS-GTR-371 for the Rothermel equations). Cite the source in a comment when adding one.
