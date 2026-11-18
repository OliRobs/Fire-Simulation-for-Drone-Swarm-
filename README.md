# Wildfire simulator

Wildfire spread on real terrain using Rothermel's surface fire spread model, driven by USGS 1 m lidar
elevation (terrain) and LANDFIRE data (Scott & Burgan 40 fuel models, canopy, operational roads). The fire runs on
10 m cells over a 20 km x 20 km area. A drone swarm surveys round the fire:
sweep drones fly lanes at right angles to the fire inside the perimeter, each in its own sector, and ACO drones
share an ant-colony-optimised tour outside it, all kept clear of the fire. The fleet grows and shrinks with the
fire so every sweep lap and every revisit on the tour stays within 20 minutes (see "Drone survey").

## Build

Needs CMake 3.16+, **SFML 2.5 – 2.6** (not SFML 3, the API changed) and GDAL 3.x.

    cmake -S . -B build && cmake --build build
    ./build/wildfire        # run from the repo root (reads the .tif files there, in ./landfire or ./dem)
    ./build/wildfire --data D:/gis/landfire --centre -117.87 34.25   # another folder / window centre
    ./build/wf_tests        # headless tests

To build just the simulation core and the tests, without SFML or GDAL:

    cmake -S . -B build -DWF_BUILD_APP=OFF && cmake --build build && ctest --test-dir build

Windows: CMake looks for SFML in `C:\SFML-2.6.1` and GDAL in `C:\gdal` (the same folders as
`UAVFireMap.vcxproj`). If yours are elsewhere, point CMake at them:

    cmake -S . -B build -DWF_SFML_ROOT=D:/libs/SFML-2.6.1 -DWF_GDAL_ROOT=D:/libs/gdal

The SFML DLLs are copied next to `wildfire.exe` after each build; GDAL's `bin` folder needs to be on `PATH`.

Controls:

| Input | Action |
|---|---|
| Left-click | Ignite a cell |
| Mouse wheel | Zoom in/out about the cursor (up to 200x, 10 cells = 100 m across) |
| Right- or middle-drag | Pan |
| `+` / `-`, arrow keys | Zoom about the centre / pan |
| `F` | Zoom to fit the fire (burning and burned cells) |
| `0` or Home | Whole map |
| `P`, `V` | Toggle the perimeter ring, drone mission overlay |
| `R` | Show/hide roads (hidden by default) |

Cell grid lines appear once cells are about 8 pixels or more on screen. The title bar shows simulated time
(one frame = 10 s), the zoom, the number of sweep drones (with the longest lap) and ACO drones (with the revisit
time), drones flying back to base, and a warning if the fleet is at its maximum and still over a time limit.
The mission overlay (`V`) tints each sweep sector in its drone's colour, draws the sweep lanes, and shows the
Voronoi regions and points of interest of the ACO tour outside the perimeter.

## Input data (USGS 1 m DEM + LANDFIRE)

Put the `.tif` files either loose in the repo root (the working folder when run from Visual Studio or as
`./build/wildfire`) or in a `landfire/` folder (sub-folders are searched), or pass `--data <folder>`. The DEM
tiles can also go in a `dem/` folder. Layers are recognised by name, so any LANDFIRE version year works:

    USGS_1M_11_x41y379_*.tif  USGS_1M_11_x41y380_*.tif  USGS_1M_11_x42y379_*.tif  USGS_1M_11_x42y380_*.tif
    LF2020_Elev_CONUS.tif   LF2020_SlpD_CONUS.tif   LF2025_FBFM40_CONUS.tif   LF2023_Roads_CONUS.tif
    LF2025_CC_CONUS.tif     LF2025_CH_CONUS.tif     LF2025_CBH_CONUS.tif      LF2025_CBD_CONUS.tif

The `.tif` files are git-ignored (LANDFIRE ~120 MB per layer, DEM tiles a few hundred MB each, too big for
GitHub). The startup log lists which layers were found.

### Terrain: USGS 3DEP 1 m DEM tiles

The four tiles x41-42 y379-380 (UTM zone 11, NAD83) cover easting 410-430 km, northing 3780-3800 km: about
34.157-34.339 N, 117.978-117.759 W, in the San Gabriel Mountains. The tile name gives its north-west corner in
10 km units (USGS tiles also carry a few metres of overlap buffer).

When tiles are found they replace LANDFIRE's 30 m elevation and slope:

- **Grid.** The tiles are joined into one mosaic (an in-memory GDAL VRT) and the simulation grid is a
  GRID_WIDTH x GRID_HEIGHT window of 10 m cells cut from it in its UTM zone (already north-up), centred on the
  tiles: 2000 x 2000 cells = exactly the four tiles. `--centre` moves the window.
- **Terrain per cell.** Slope is worked out on every 1 m pixel (Horn's method) and averaged into the cell, so a
  cell crossing a gully or steep bank keeps that steepness (slope from 10 m averaged elevations would smooth it
  away, and Rothermel's slope factor grows with tan^2 of the slope). Elevation is the mean of the 100 pixels; the
  upslope direction (aspect) is that of the mean gradient. See `aggregateFineTopography` in `topography.cpp`.
- **Memory.** The 400 million pixels are read in strips of 50 cell rows (~40 MB each), never all at once.
  Loading takes roughly 20 s.
- **Why not 1 m cells?** 400 million cells would need ~16 GB and tiny time steps, and 1 m is below the scale the
  Rothermel model is meant for (a uniform fuel bed); the LANDFIRE fuel data is 30 m anyway.
- Cells the tiles don't cover (only if the window is moved off them) fall back to the LANDFIRE elevation and slope.

The LANDFIRE layers are warped onto the same 10 m cells: fuel by most common value, roads by max, canopy by
average. They are still 30 m data, so each 30 m pixel fills about 3 x 3 cells. The LANDFIRE download must cover
the tiles (check its GeoJSON footprint).

| Layer | Used for | Notes |
|---|---|---|
| 1 m DEM tiles (`USGS_1M_`, m) | Simulation grid, elevation, slope and aspect | See above. Takes priority over the LANDFIRE elevation and slope. |
| Elevation (`_Elev_`, m) | Without DEM tiles: simulation grid, background shading, **aspect** (which way is uphill) | The download is 30 m pixels over ~264 x 207 km, so a window of GRID_WIDTH x GRID_HEIGHT cells (2000 x 10 m = 20 km) is cut out, centred on the raster centre or on `--centre <lon> <lat>`. Not read when the DEM tiles cover the whole window. |
| Slope (`_SlpD_`, degrees) | Without DEM tiles: Rothermel slope factor | Missing cells (or a missing layer) fall back to slope computed from the elevation. |
| FBFM40 (`FBFM40`) | Fuel model per cell | Scott & Burgan model numbers, looked up in `fuel_models.cpp` (Table 7 of RMRS-GTR-153). 91-99 (urban, water, agriculture, bare...) and NoData don't burn. |
| Roads (`_Roads_`) | Fuel breaks (not drawn unless `R` is pressed) | Classes 20 primary, 21 secondary, 22 tertiary, 23 thinned. Road cells don't burn (`ROADS_ARE_FUEL_BREAKS`); diagonal road lines hold too. |
| Canopy cover (`_CC_`, %) | Wind sheltering | With canopy height and base height gives the crown fill portion, which sets the sheltered wind adjustment factor (Andrews 2012). |
| Canopy height (`_CH_`, m x 10) | Wind sheltering | |
| Canopy base height (`_CBH_`, m x 10) | Crown fire initiation, crown ratio | Van Wagner (1977) torching threshold. |
| Canopy bulk density (`_CBD_`, kg/m3 x 100) | Active crown fire | Crown fire is active only if its spread rate exceeds 3.0 / CBD m/min. |

Canopy values are stored by LANDFIRE as scaled integers (0 = no forest) and converted on load (RMRS-GTR-220).

**North-up grid.** LANDFIRE's CONUS Albers projection (EPSG:5070) only points north at 96 W; elsewhere its grid
is rotated (about 13 degrees in Southern California), which made the map, roads and wind direction tilted. So
the window is laid out north-up in the local UTM zone (NAD83, picked from the window centre) and every layer is
warped onto it: average for elevation, slope and canopy, most common value for fuel, max for roads (so a
one-pixel road stays an unbroken fuel break). Set `ALIGN_NORTH_UP = false` in `config.h` to use the Albers
pixel grid unchanged instead.
Every layer is warped onto the elevation layer's grid (average for elevation and slope, most common value for
fuel and roads), so layers with different extents still line up. Without a LANDFIRE folder the program falls
back to the bundled UK `elevation_map.tif` (used whole, as before) with fuel model GR2 everywhere.

## Fire spread: Rothermel surface fire model

`rothermel.cpp` implements Rothermel (1972) with Albini's (1976) revisions, as set out in Andrews (2018),
RMRS-GTR-371 (Table 6), for the multi-class Scott & Burgan fuel models, including the dynamic transfer of cured
herbaceous load from live to dead:

    R = I_R xi (1 + phi_w + phi_s) / (rho_b epsilon Q_ig)

- Per fuel model: reaction intensity, propagating flux ratio, heat sink and the no-wind no-slope rate R_0.
- Per cell: wind and slope are added as vectors (wind in any direction relative to the slope) to give the head
  fire rate and direction, and an elliptical fire shape from the effective wind speed (RMRS-GTR-371 section 6,
  as in BehavePlus). The ellipse gives the spread rate towards each of the 8 neighbours.
- On the grid: the travel time between two cells uses half a step at each cell's rate, and fire arrival times are
  propagated in time order, so the front moves at the Rothermel rates whatever the frame rate.
- A cell shows as burning while the front crosses it plus the flame residence time t_r = 384/sigma.

### Canopy and crown fire (`crown_fire.cpp`)

Under a forest canopy the midflame wind is reduced with the sheltered wind adjustment factor (Andrews 2012,
crown fill portion from CC, CH and CBH). Crown fire follows the Finney (1998) method that FARSITE and FlamMap
use by default:

1. **Torching** (Van Wagner 1977): if the surface fireline intensity (Byram, I_B = I_R t_r R) reaches
   I_0 = (0.010 CBH (460 + 25.9 FMC))^1.5 kW/m, the canopy ignites.
2. **Crown spread rate** (Rothermel 1991): 3.34 x the surface rate of original fuel model 10 with midflame wind
   0.4 x the 20-ft wind.
3. **Active vs passive**: if that rate is at least 3.0 / CBD m/min the cell is an active crown fire and spreads at
   the crown rate (or the surface rate if that is faster); otherwise it only torches and spreads at the surface
   rate.

Burning cells are coloured by fire type: orange-red surface fire, orange torching, yellow active crown fire.
The startup log reports how many cells would torch or crown under the current weather.
Crown fire spread is not affected by the surface WAF: the crown is in the open wind.

Check: with Scott & Burgan's standard conditions (D2L2 moisture, 5 mi/h midflame wind, flat) all 40 fuel models
land in the spread-rate class the paper gives them, and packing ratios and characteristic SAV match the
published values (`wf_tests`).

### Data the model still needs (not in the LANDFIRE layers)

These are fixed placeholder values in `config.h` / `FireWeather` for now:

| Input | Current value | Where to get it |
|---|---|---|
| **Dead fuel moisture** 1-h, 10-h, 100-h | 6 / 7 / 8 % (Scott & Burgan scenario D2) | Nearest RAWS weather station / NFDRS output, or computed from temperature and relative humidity (e.g. Nelson dead fuel moisture model). The single biggest driver of spread rate. |
| **Live fuel moisture** herbaceous, woody | 60 / 90 % (scenario L2, two-thirds cured) | National Fuel Moisture Database samples, or a curing / greenness estimate (RMRS-GTR-153 Table 2 converts % cured to herb moisture). Controls the dynamic grass models strongly. |
| **Foliar moisture** of the canopy | 100 % | Seasonal values from the literature (normal range 75-150 %, RMRS-GTR-220). Lower values make torching much easier. |
| **Wind speed** at 20 ft | 20 km/h | RAWS observations or a forecast (NWS / HRRR). Decides whether crown fire becomes active. |
| **Wind direction** | from 225 deg (SW) | As above. Real wind is bent by terrain: a gridded wind (e.g. WindNinja run on the LANDFIRE DEM) would replace the single value. |

Aspect is not needed as a separate layer: it is computed from the elevation (Horn's method, `topography.cpp`),
from the 1 m DEM when the tiles are present.

## Layout (src/)

| Module | Responsibility |
|---|---|
| config.h, types.h | constants, neighbour tables, Vec2, Cell |
| raster_io, terrain | GDAL loading: finds the DEM tiles and LANDFIRE layers, cuts the window, reads the 1 m DEM in strips, warps every layer onto the grid (only code using GDAL) |
| viewport | Zoom and pan maths (headless, tested); the renderer turns it into an sf::View |
| fuel_models | The 40 Scott & Burgan fuel models (Table 7 of RMRS-GTR-153), looked up by FBFM40 code |
| rothermel | Rothermel surface spread: fuel bed, wind and slope factors, wind adjustment (open and under canopy), head fire direction and ellipse |
| crown_fire | Van Wagner torching threshold, active crown fire threshold, fireline intensity and flame length |
| topography | slope and aspect (upslope direction) from the DEM; 10 m cell terrain from 1 m pixels |
| fire_model | FireSimulation: per-cell Rothermel rates, fire arrival times in time order |
| distance_transform | exact Euclidean distance transform; optionally range-limited to a window round the fire so the per-frame cost follows the fire size |
| safe_planner | fire-safe legs: line of sight, A*, string pulling; reports "no safe path" instead of guessing |
| perimeter | committed fire + sectional perimeter updates; traces the perimeter line into ordered loops |
| voronoi, boustrophedon, aco | the three planning algorithms (boustrophedon: sweep lanes at right angles to the fire, gap fill, ownership map) |
| mission_planner | zones, Voronoi regions (ACO points of interest), sweep lanes and sectors, fleet sizing, sector-to-drone matching, ACO tour |
| drone, mission_controller | drone movement and safety; re-planning, re-routing blocked drones, route assignment; adds and retires drones |
| renderer | all SFML drawing |
| main.cpp | window, events, frame loop |

Library layering: wf_core (everything except GDAL/SFML) <- wf_io (GDAL), wf_render (SFML) <- wildfire.

## Key settings (config.h)

The grid is GRID_WIDTH x GRID_HEIGHT = 2000 x 2000 cells of SIM_CELL_SIZE_M = 10 m. The drone and perimeter
distances are set in metres and converted to cells with `cellsFromMetres`, so they keep the same ground
distances they had on the original 30 m cells:

| Setting | Value | Meaning |
|---|---|---|
| PERIMETER_DISTANCE | 1200 m (120 cells) | Perimeter ring distance from the committed fire |
| PERIMETER_UPDATE_DISTANCE | 600 m (60 cells) | A perimeter section moves once the fire gets this close to it |
| MIN_FIRE_DISTANCE | 300 m (30 cells) | Survey standoff for sweeps and points of interest |
| HARD_MIN_FIRE_DISTANCE | 90 m (9 cells) | Flight legs may cut inside the standoff, but never closer than this |
| SEARCH_DISTANCE | 2100 m (210 cells) | Outer edge of the ACO search area |
| SENSOR_RADIUS | 60 m (6 cells) | Sensor swath half-width |
| DRONE_SPEED_MPS | 16.5 m/s (16.5 cells per 10 s tick) | Drone ground speed |
| TARGET_REGION_CELLS | 630 x 630 m (4050 cells) | Voronoi region size target |
| LANE_SPACING | 120 m (12 cells) | Most distance between sweep lanes along the perimeter line |
| MAX_SWEEP_LAP_S | 20 min | Longest sweep lap (closed route of one sector, incl. the leg back to the start) |
| SWEEP_LAP_MARGIN | 0.9 | Sectors are first sized to 90 % of the lap limit, leaving room for fire detours |
| MIN_SWEEP_DRONES, MAX_SWEEP_DRONES | 3, 10 | Range of the sweep fleet |
| ACO_REVISIT_S | 20 min | Longest wait between visits to any point on the ACO tour |
| MIN_ACO_DRONES, MAX_ACO_DRONES | 1, 10 | Range of the ACO fleet |
| SHRINK_MARGIN, SHRINK_REPLANS | 0.85, 3 | A drone is removed only after 3 re-plans in a row on which one fewer would stay within 85 % of the limit |
| SWEEP_DRONES, ACO_DRONES | 3, 5 | Fleet at the start; it then follows the fire |

## Drone survey

**Zones.** Cells under 1.2 km from the committed fire and at least 300 m from the live fire form the inner band,
swept by the sweep drones. From 1.2 km out to 2.1 km is the outer zone, covered by the ACO tour through one point
of interest per Voronoi region.

**Sweep lanes.** The 1.2 km line is traced from the distance field into ordered closed loops (marching squares;
one loop per separate fire, all anticlockwise round the fire). Lanes are seeded evenly along each loop, at most
120 m apart, and run straight in down the distance gradient, i.e. at right angles to the nearest fire edge,
keeping the parts inside the band (normally down to the 300 m standoff). On convex stretches lanes converge as
they go in; that overlap is accepted. In concave bays they spread apart, so short gap-fill lanes are added
until every inner cell is within the sensor radius of a lane.

**Sectors.** Each sweep drone owns a contiguous run of lanes along the perimeter line, flown as a closed loop,
lanes alternately inward and outward. Borders fall halfway between two lanes, so each sector's first lane is
60 m inside it. Drones may fly across other drones' sectors in transit, on re-routes, on fire detours and on the
way home; only survey lanes are exclusive. After each re-plan sectors are matched to drones by the land they
shared before, so a drone keeps its side of the fire, and each drone joins its route at the nearest waypoint.

**Sweep fleet size.** Lap time is the length of the closed fire-safe route divided by the drone speed, which is
also the longest wait between visits to any lane. At each re-plan the planner finds the fewest sectors whose
straight-line laps fit in 90 % of 20 minutes, splits the lanes so the longest lap is as short as possible, then
builds and times the real routes and adds a drone if a lap is still over 20 minutes. Round fires need about:

| Fire across | Perimeter line | Sweep drones |
|---|---|---|
| New fire | 7.6 km | 4 |
| 1 km | 10.7 km | 6 |
| 2 km | 13.8 km | 8 |
| 3 km | 17.0 km | 10 (maximum) |

A long, wind-driven fire has a longer perimeter line than a round one of the same width and needs more.

**ACO tour.** The tour is always flown anticlockwise. After a re-plan each ACO drone locks back on at the
nearest node of the new tour (or the next one if it has just passed it) and carries on the same way, so drones
never turn back. New drones join in the middle of the largest gaps; when the fleet shrinks, the drone in the most
crowded spot goes home, keeping the rest as evenly spaced as possible.

**ACO fleet size.** The ACO drones share one tour, evenly spaced, each flying all of it. Their count is the time
one drone takes to fly the tour divided by 20 minutes, rounded up, so every point is passed at least every 20
minutes. One drone covers a new fire's tour (about 17 km); a 2 km fire needs two.

**Growing and shrinking.** New drones launch from the launch point as soon as they are needed. A drone is removed
only when one fewer would keep every lap (or the revisit time) within 85 % of the limit on 3 re-plans in a row,
at least 30 simulated minutes; it then flies a fire-safe route back to the launch point and leaves the fleet
there (or holds if no safe route exists). Which drone retires is chosen in
`MissionController::chooseDroneToRetire`, ready for battery level to be taken into account. At the maximum of
10 drones per role, laps may go over the limit and the title bar shows a warning. Collision avoidance between
drones, batteries and recharging are not modelled yet.

## Drone safety

- Every planned leg is checked against the fire; if a waypoint can only be reached by crossing the fire it is
  dropped from the route rather than flown in a straight line.
- Between re-plans the fire keeps growing, so each segment is re-checked as the drone reaches it. A blocked
  drone is re-routed straight away, or hovers if nothing is reachable, retrying every 15 frames.
- `keepClearOfFire` is the last-resort push away from the fire.
