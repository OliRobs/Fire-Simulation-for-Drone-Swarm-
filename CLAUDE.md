# Wildfire simulator with drone survey

C++17 cellular-automaton wildfire on real terrain (GDAL GeoTIFFs, drawn with SFML), plus a drone swarm that
surveys round the fire. README.md has the module table and input-data notes.

## Build and test

- Full build: `cmake -S . -B build && cmake --build build`, then run `./build/wildfire` from the repo root
  (it loads `elevation_map.tif` and optional `fuel_map.tif` from the working directory).
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
- `elevation_map.tif` is British National Grid (EPSG:27700), 5 km square. The DEM defines the simulation
  grid; every other layer is warped onto it. LANDFIRE data covers only the US, so it doesn't overlap this DEM.
- Fuel rasters hold FBFM40 category codes: look them up through `fuel_models.cpp`, never average or clamp them.

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
- Physics and fuel constants come from the papers in the claude.ai project (Scott & Burgan 2005; Rothermel;
  Alexandridis et al. 2008 for slope). Cite the source in a comment when adding one.
