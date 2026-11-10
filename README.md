# Wildfire simulator

Wildfire spread on real terrain using Rothermel's surface fire spread model, driven by LANDFIRE data
(elevation, slope, Scott & Burgan 40 fuel models, operational roads). A drone swarm surveys round the fire:
boustrophedon sweeps inside the perimeter and an ant-colony-optimised tour outside it, both kept clear of
the fire.

## Build

Needs CMake 3.16+, **SFML 2.5 – 2.6** (not SFML 3, the API changed) and GDAL 3.x.

    cmake -S . -B build && cmake --build build
    ./build/wildfire        # run from the repo root (reads the LANDFIRE .tif files there or in ./landfire)
    ./build/wildfire --data D:/gis/landfire --centre -118.07 34.29   # another folder / window centre
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
| Mouse wheel | Zoom in/out about the cursor (up to 50x, 10 cells across) |
| Right- or middle-drag | Pan |
| `+` / `-`, arrow keys | Zoom about the centre / pan |
| `F` | Zoom to fit the fire (burning and burned cells) |
| `0` or Home | Whole map |
| `P`, `V` | Toggle the perimeter ring, drone mission overlay |
| `R` | Show/hide roads (hidden by default) |

Cell grid lines appear once cells are about 8 pixels or more on screen. The title bar shows simulated time
(one frame = 10 s) and the zoom.

## Input data (LANDFIRE)

Put the LANDFIRE `.tif` files either loose in the repo root (the working folder when run from Visual Studio
or as `./build/wildfire`) or in a `landfire/` folder (sub-folders are searched), or pass `--data <folder>`.
Layers are recognised by name, so any LANDFIRE version year works:

    LF2020_Elev_CONUS.tif   LF2020_SlpD_CONUS.tif   LF2025_FBFM40_CONUS.tif   LF2023_Roads_CONUS.tif
    LF2025_CC_CONUS.tif     LF2025_CH_CONUS.tif     LF2025_CBH_CONUS.tif      LF2025_CBD_CONUS.tif

The `.tif` files are git-ignored (about 120 MB each, too big for GitHub). The startup log lists which layers
were found.

| Layer | Used for | Notes |
|---|---|---|
| Elevation (`_Elev_`, m) | Simulation grid, background shading, **aspect** (which way is uphill) | The download is 30 m pixels over ~264 x 207 km, so a window of GRID_WIDTH x GRID_HEIGHT cells (500 x 30 m = 15 km) is cut out, centred on the raster centre or on `--centre <lon> <lat>`. |
| Slope (`_SlpD_`, degrees) | Rothermel slope factor | Missing cells (or a missing layer) fall back to slope computed from the elevation. |
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

Aspect is not needed as a separate layer: it is computed from the elevation (Horn's method, `topography.cpp`).

## Layout (src/)

| Module | Responsibility |
|---|---|
| config.h, types.h | constants, neighbour tables, Vec2, Cell |
| raster_io, terrain | GDAL loading: finds the LANDFIRE layers, cuts a north-up UTM window, warps every layer onto it (only code using GDAL) |
| viewport | Zoom and pan maths (headless, tested); the renderer turns it into an sf::View |
| fuel_models | The 40 Scott & Burgan fuel models (Table 7 of RMRS-GTR-153), looked up by FBFM40 code |
| rothermel | Rothermel surface spread: fuel bed, wind and slope factors, wind adjustment (open and under canopy), head fire direction and ellipse |
| crown_fire | Van Wagner torching threshold, active crown fire threshold, fireline intensity and flame length |
| topography | slope and aspect (upslope direction) from the DEM |
| fire_model | FireSimulation: per-cell Rothermel rates, fire arrival times in time order |
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
