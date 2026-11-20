# Wildfire simulator with drone survey

C++17 cellular-automaton wildfire on real terrain (GDAL GeoTIFFs, drawn with SFML), plus a drone swarm that
surveys round the fire. README.md has the module table and input-data notes.

## Build and test

- Full build: `cmake -S . -B build && cmake --build build`, then run `./build/wildfire` from the repo root
  (it loads the LANDFIRE layers from the repo root or `landfire/`, falling back to `elevation_map.tif`).
- Tests: `./build/wf_tests` (or `ctest --test-dir build`). They are headless and take about 25 s (most of
  it the end-to-end mission test on the 2000 x 2000 grid).
- Core and tests only, without SFML or GDAL: add `-DWF_BUILD_APP=OFF`.
- The owner builds day to day with Visual Studio: `UAVFireMap.vcxproj`, SFML in `C:\SFML-2.6.1`, GDAL in
  `C:\gdal`. The project is x64 only (the SFML and GDAL builds are 64-bit). Debugger settings (working
  directory, PATH to the DLLs) live in the .vcxproj itself, because `.vcxproj.user` is not tracked.
  Adding or removing a source file means updating `CMakeLists.txt`, `UAVFireMap.vcxproj` AND
  `UAVFireMap.vcxproj.filters`.
- GDAL on Windows is the GISInternals build: DLLs in `C:\gdal\bin`, PROJ's `proj.db` in `C:\gdal\bin\proj9\share`,
  GDAL data in `C:\gdal\bin\gdal-data`. `registerRasterDrivers` points PROJ and GDAL there when `PROJ_DATA` /
  `GDAL_DATA` aren't set; without proj.db no LANDFIRE layer can be reprojected onto the DEM grid.
- CI (`.github/workflows/build.yml`): full build and tests on Ubuntu 24.04, plus an MSVC build of the core and
  tests on Windows. Both must pass.

## Constraints

- SFML stays on 2.6. The code uses the SFML 2 API (`sf::Uint8`, `VideoMode(w, h)`, `pollEvent(event)`);
  don't port to SFML 3, which would break the Visual Studio setup. vcpkg's `sfml` port is SFML 3, so don't
  use it.
- Library layering: only `raster_io` and `terrain` use GDAL; only `renderer` and `main.cpp` use SFML.
  Everything else is in `wf_core`, which must build without either so the tests stay headless.
- The values in `config.h` (perimeter distances, fire standoffs, drone counts, speeds) are design decisions.
  Ask before changing them. Distances and speeds are set in metres and converted with `cellsFromMetres`, so
  changing SIM_CELL_SIZE_M keeps their ground values.
- Terrain comes from USGS 3DEP 1 m DEM tiles (`USGS_1M_11_x41-42y379-380_*.tif`, UTM zone 11, San Gabriel
  Mountains). They define the simulation grid: GRID_WIDTH x GRID_HEIGHT = 2000 x 2000 cells of SIM_CELL_SIZE_M =
  10 m, i.e. exactly the four tiles. Each cell's elevation, slope (mean of the 1 m slopes) and aspect come from the
  1 m pixels inside it (`aggregateFineTopography`), read in strips so the 1 m data is never all in memory. Do not
  run the fire on 1 m cells (400 M cells, ~16 GB, and below Rothermel's scale).
- Other input data is LANDFIRE (US only, EPSG:5070 Albers, 30 m): FBFM40 fuel models, operational roads and the
  canopy layers CC (%), CH (m x 10), CBH (m x 10), CBD (kg/m3 x 100), plus elevation and slope (degrees), which are
  used only where there are no DEM tiles. Files are found by name in `--data`, else `landfire/` (DEM tiles also in
  `dem/`), else loose in the working folder (the repo root). Without DEM tiles the LANDFIRE elevation layer defines
  the grid, laid out north-up in the local UTM zone (`ALIGN_NORTH_UP`; the Albers grid is rotated ~13 deg in
  California). Every LANDFIRE layer is warped onto the grid (fuel Mode, roads Max).
  The .tif files are git-ignored. `elevation_map.tif` (UK, EPSG:27700, 5 km) is only the no-data fallback.
- Fuel rasters hold FBFM40 category codes: look them up through `fuel_models.cpp`, never average or clamp them
  (resample with Mode; roads with Max so thin roads stay continuous).
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
  `SEARCH_DISTANCE`. Only the ACO region (the outer zone) is split into Voronoi regions, which give the ACO
  points of interest; the inner band is not segmented with Voronoi cells.
- Sweep lanes run straight in from the perimeter line (traced from the committed-fire distance field), at right
  angles to the nearest fire edge, seeded at most `LANE_SPACING` (120 m) apart. Overlap where lanes converge is
  accepted; gap-fill lanes cover the cells lanes miss in concave bays.
- Each sweep drone owns one sector: a contiguous run of lanes along the perimeter line, borders halfway between
  lanes. Only survey lanes are exclusive; drones may cross other sectors in transit, re-routes, fire detours and
  on the way home. Sectors are matched to drones by the land they shared before (stable ownership); after a
  re-plan a drone joins its route at the nearest waypoint and carries on in route order.
- Fleet size follows the fire: sweep laps (closed route incl. the closing leg) within `MAX_SWEEP_LAP_S`; ACO count
  = ceil(one-drone tour time / `ACO_REVISIT_S`). 3-10 sweep and 1-10 ACO drones. Growing is immediate (launch from
  the launch point); a drone is retired only after `SHRINK_REPLANS` re-plans in a row within `SHRINK_MARGIN` of the
  limit with one fewer, and flies a fire-safe route home before it is removed. Drones are found by role and id,
  never by position in the list. Collision avoidance and batteries are out of scope for now.
- The ACO routing only changes when the perimeter does: while no perimeter section has moved, a re-plan keeps the
  Voronoi regions, points of interest and tour order (`SweepFleet::perimeterUnchanged`, from
  `DistanceTransform::version()` of the perimeter's field) and only rebuilds the fire-safe legs.
- ACO drones first fly to the closest point on the tour, then follow it, evenly spaced and in one direction:
  the tour is always handed out anticlockwise on screen. After a re-plan, a drone already on the tour locks back
  on at the nearest node (the next one if it has just passed it) and carries on anticlockwise; it is never sent
  back or across to another slot. New drones fill the largest gaps; when the fleet shrinks, the drone in the most
  crowded spot retires. Even spacing is restored by slowing down, never by turning back: an ACO drone closer than
  an even share of the tour behind the drone ahead flies slower (down to `ACO_SPACING_MIN_SPEED`) until the gap
  has opened up.
- The perimeter updates in sections, only where the fire front has come within `PERIMETER_UPDATE_DISTANCE`.

## Conventions

- Comment all code clearly: a header block above each function saying what it does and why, and comments on
  non-obvious steps. Match the existing style (`// ----` section banners).
- Any change to fire spread, planning or drone movement needs a test in `tests/tests.cpp`. Run `wf_tests`
  before committing.
- Work on a branch and open a pull request; don't push to `main`.
- Physics and fuel constants come from the papers in the claude.ai project (Scott & Burgan 2005 for the fuel
  models; Andrews 2018 RMRS-GTR-371 for the Rothermel equations). Cite the source in a comment when adding one.
