# Wildfire simulator

Wildfire spread on real terrain using Rothermel's surface fire spread model, driven by USGS 1 m lidar
elevation (terrain) and LANDFIRE data (Scott & Burgan 40 fuel models, canopy, operational roads). The wind is
adjusted over the terrain by a mass-conserving wind solver (WindNinja's method), so each cell has its own wind
speed and direction (see "Terrain wind"). The fire runs on 10 m cells over a 12 km x 10 km area
(1200 x 1000 cells). A drone swarm surveys round the fire:
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

**Benchmark.** `./build/wf_bench [frames]` (default 2400 frames, 200 simulated minutes) runs the frame loop
headless on a flat grass map lit at the centre and prints the time of each phase per frame, the re-plan average
and worst case, and a checksum of the final state (fire, distance fields, plan, routes, drones). CI runs it with
`--expect` and fails if the checksum differs from `tools/bench_checksum.txt`: an optimisation must leave it
unchanged; a change meant to alter results (including a `config.h` design value) updates that file in the same
pull request.

`./build/wf_bench --wind` times the terrain wind solver on its own, at the app's size (10 m cells, a 120 x 100 x 20
wind mesh) on a synthetic ridge-and-valley map, and prints a checksum of every cell's solved wind. CI checks it
against `tools/wind_checksum.txt` in the same way.

Windows: CMake looks for SFML in `C:\SFML-2.6.1` and GDAL in `C:\gdal` (the same folders as
`UAVFireMap.vcxproj`). If yours are elsewhere, point CMake at them:

    cmake -S . -B build -DWF_SFML_ROOT=D:/libs/SFML-2.6.1 -DWF_GDAL_ROOT=D:/libs/gdal

The SFML DLLs are copied next to `wildfire.exe` after each build; GDAL's `bin` folder needs to be on `PATH`.
GDAL also needs PROJ's database, `proj.db`, to reproject the LANDFIRE layers onto the DEM grid. If the
`PROJ_DATA` environment variable isn't set, the program looks in `C:\gdal\bin\proj9\share` (the GISInternals
layout) and a few other usual places, and says at startup where it found it. Without it you get
"Cannot find proj.db" and every LANDFIRE layer is skipped (fuel model GR2 everywhere, no roads or canopy).

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
| `W` | Show/hide the terrain wind arrows (hidden by default; see "Terrain wind") |

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
  tiles. GRID_WIDTH x GRID_HEIGHT is 1200 x 1000 cells (12 km x 10 km), so the window covers the middle of
  the 20 km x 20 km mosaic: about easting 414-426 km and northing 3785-3795 km, leaving 4 km of tile east and
  west and 5 km north and south. `--centre` moves the window; it is still cut on the DEM's own grid (edges on
  whole metres), and cells moved off the tiles use the LANDFIRE elevation and slope.
- **Terrain per cell.** Slope is worked out on every 1 m pixel (Horn's method) and averaged into the cell, so a
  cell crossing a gully or steep bank keeps that steepness (slope from 10 m averaged elevations would smooth it
  away, and Rothermel's slope factor grows with tan^2 of the slope). Elevation is the mean of the 100 pixels; the
  upslope direction (aspect) is that of the mean gradient. See `aggregateFineTopography` in `topography.cpp`.
- **Memory.** The window's 120 million 1 m pixels are read in strips of 50 cell rows (~24 MB each), never all
  at once. Loading takes several seconds (roughly 20 s when the window was the whole 20 km mosaic).
- **Why not 1 m cells?** 120 million cells for this window would need ~5 GB and tiny time steps, and 1 m is below the scale the
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
| **Wind direction** | from 225 deg (SW) | As above. This is the open wind above the terrain; the terrain wind solver turns it along valleys and round hills cell by cell (see "Terrain wind"). |

Aspect is not needed as a separate layer: it is computed from the elevation (Horn's method, `topography.cpp`),
from the 1 m DEM when the tiles are present.

## Terrain wind (`wind_solver.cpp`)

The weather gives one open wind: 20 km/h from 225 deg, 20 ft (6.1 m) above the vegetation (`config.h`). Real
wind does not stay the same across a mountain range. It speeds up over ridges and summits, slows in hollows and
on the slopes facing into it, and is turned along valleys and round hills. Fire behaviour follows the wind, so
`solveTerrainWind` works out a separate 20-ft wind for every cell before the fire starts. It uses the method of
WindNinja's default "conservation of mass" solver (Forthofer et al. 2014; github.com/firelab/windninja, analysed
in the project doc "How WindNinja Calculates Wind"). The idea: take a first guess at the wind that ignores the
terrain, then make the smallest change to it that stops air piling up against hillsides or vanishing out of
valleys, i.e. that conserves mass.

### How it works

1. **Wind mesh.** The fire grid is grouped into 100 m columns (10 x 10 fire cells each, `WIND_COLUMN_SIZE_M`).
   Each column gets a stack of 20 nodes (`WIND_LAYERS`) from its mean ground elevation up to a flat top 2.5 km
   above the highest summit (`WIND_TOP_ABOVE_SUMMIT_M`). The layers are thin near the ground (a few metres) and
   each is 1.3 x thicker than the one below, so the bottom layers follow the terrain and the top is flat. On the
   default 1200 x 1000 grid that is 120 x 100 x 20 = 240,000 nodes. Solving on the 10 m fire cells instead would
   mean ~24 million unknowns, and ridge and valley effects are hundreds of metres across anyway.
2. **First guess.** The open wind is spread up every column with a neutral logarithmic profile,

       U(z) = U20 * ln((z - d) / z0) / ln((6.1 m + h - d) / z0)

   where h is the vegetation height (the LANDFIRE canopy height where there are trees, else the fuel bed depth of
   the cell's fuel model), the roughness length is z0 = 0.13 h (at least 1 cm) and the displacement height is
   d = 0.64 h. These are the same ratios the Albini & Baughman wind adjustment factor in `rothermel.cpp` is built
   on, so the wind field and the midflame wind agree. Below d + 7 z0 the speed falls linearly to zero; above the
   boundary-layer top, 0.2 u* / f, it stops increasing. The first guess blows horizontally everywhere, at the
   same height above the ground, so it blows straight into hillsides: it has divergence (mass errors).
3. **The adjustment equation.** Find the wind (u, v, w) closest to the first guess (u0, v0, w0) with
   div(u) = 0. Minimising (u - u0)^2 + (v - v0)^2 + (w - w0)^2 / alpha^2 over the domain, with a Lagrange
   multiplier Phi for the constraint, gives a single linear equation (Sasaki 1958, Sherman 1978):

       d/dx(Rx dPhi/dx) + d/dy(Ry dPhi/dy) + d/dz(Rz dPhi/dz) = -(du0/dx + dv0/dy + dw0/dz)
       Rx = Ry = 1/2,   Rz = alpha^2 / 2

   and the corrected wind is u = u0 + Rx dPhi/dx, v = v0 + Ry dPhi/dy, w = w0 + Rz dPhi/dz. The right-hand side
   is the first guess's mass error; Phi is the "pressure" that pushes the air round and over the terrain to
   cancel it. alpha (`WIND_STABILITY_ALPHA`) sets how easily air moves vertically compared with horizontally.
4. **Finite elements.** Each mesh cell is an 8-node (trilinear) hexahedron, integrated with 2 x 2 x 2 Gauss
   points. Because the mesh is structured, each node is coupled to at most its 27 neighbours, so the matrix is
   stored as 27 numbers per node. Boundary conditions follow WindNinja:
   - Ground: nothing imposed (the natural condition), which means no air flows through the ground.
   - Four sides and the top: Phi = 0, so air flows in and out freely and the first guess is kept there. The
     edge cells of the map therefore keep the open wind; the effect of the terrain builds up within a column or
     two of the edge.

   The matrix is symmetric and positive definite, so it is solved with conjugate gradients, preconditioned with
   symmetric Gauss-Seidel (SSOR with omega = 1, as WindNinja). It stops when the residual is 10^-4 of the
   right-hand side: about 90 iterations, 1.5 s for the full grid.
5. **Wind at the nodes.** The gradient of Phi jumps between elements, so each node takes the
   inverse-distance-weighted average of the gradients at the surrounding Gauss points (WindNinja's "stress
   smoothing"), and step 3's formulas give u, v, w.
6. **20-ft wind per fire cell.** In each column the horizontal wind is read 6.1 m above the vegetation:
   - in the lowest mesh layer, by scaling the first node's wind down the log profile;
   - higher up, by interpolating wind / ln-profile linearly between the two nodes and multiplying back.
     Interpolating along the profile rather than in a straight line means an undisturbed log profile comes back
     exactly; a straight line reads about 1 % low.

   The column winds are then interpolated bilinearly (by component) to every 10 m fire cell, giving a
   `WindField`: a speed in km/h and a unit direction per cell.

### How the fire uses it

`FireSimulation` takes the `WindField` as an optional last argument (`main.cpp` passes it when
`USE_TERRAIN_WIND` is on). For each cell it replaces the single open wind with the cell's own 20-ft wind:

- **Speed:** 20-ft wind x wind adjustment factor (fuel depth or canopy sheltering, unchanged) = midflame wind,
  which sets Rothermel's wind factor phi_w.
- **Direction:** the wind vector added to the slope vector to give the head fire direction and the spread
  ellipse, as before.
- **Crown fire:** the crown spread rate uses 0.4 x the cell's own 20-ft wind.

Without a `WindField` (`USE_TERRAIN_WIND = false`, the tests and the benchmark) every cell gets the single open
wind as before, with identical results. The startup log reports the mesh, the CG iterations and residual, the
time taken, the range of 20-ft wind speeds across the map and the mass error (RMS divergence) before and after.

### Seeing it: wind arrows (`W`)

Press `W` to draw the wind field over the map (`wind_overlay.cpp` picks the arrows, `renderer.cpp` draws them):

- **Direction:** each arrow points where the 20-ft wind blows, in the cell under its centre.
- **Length:** in proportion to the speed; the fastest wind on the map gets the longest arrow (its speed is
  printed to the console when you turn the arrows on).
- **Colour:** the speed against the open wind from the weather. White = the same as the open 20 km/h,
  shading to blue where the terrain slows it (fully blue at 30 % or more slower) and to magenta where it
  speeds it up (fully magenta at 30 % or more faster). Valleys and lee slopes come out blue, ridges magenta.
- **Spacing:** the arrows sit on a lattice fixed to the map, about 48-96 screen pixels apart at any zoom, so
  zooming in shows more of them and panning doesn't move them about. Close in, neighbouring arrows look alike
  because the wind is solved on 100 m columns.

The key does nothing (and says so in the console) when `USE_TERRAIN_WIND` is off.

### Settings (`config.h`)

| Setting | Value | Meaning |
|---|---|---|
| USE_TERRAIN_WIND | true | Solve the terrain wind; false = one open wind everywhere |
| WIND_COLUMN_SIZE_M | 100 m | Horizontal spacing of the wind mesh |
| WIND_LAYERS | 20 | Mesh nodes up each column (WindNinja's default) |
| WIND_TOP_ABOVE_SUMMIT_M | 2500 m | Height of the mesh top above the highest summit |
| WIND_STABILITY_ALPHA | 1.0 | Atmospheric stability: 1 neutral, < 1 stable (down to 0.2, clear calm night), > 1 unstable (up to 5, strong sun, light wind) |

**Why the top is higher than WindNinja's.** WindNinja sets its domain top from the cell size (column width /
400, grown over the 20 layers), which is only ~160 m above the highest summit at 100 m columns. The top is an
open boundary, so with it that low the air escapes upwards instead of speeding up over the ridges. In tests on
a 400 m high hill this cut the crest speed-up from ~20 % to ~2 %. The answer stops changing once the top is
about 2.5 km up, so that is the default; WindNinja's height is still used if it is higher.

**Stability.** With alpha < 1 (stable air, e.g. at night) vertical motion is expensive, so more air is turned
round hills and the horizontal speed changes are larger. With alpha > 1 (unstable, sunny afternoon) air goes
over the terrain and the wind changes less. On the test hill the wind on the flank is turned 7 deg round it at
alpha = 0.3, 2 deg at alpha = 1 and 0.3 deg at alpha = 3. WindNinja picks alpha from the Pasquill stability
class (sun strength, cloud and wind speed); here it is a fixed setting until weather data is wired in.

### Checks (`wf_tests`)

On 30 m cells with 600 m columns (fast; the full 10 m / 100 m resolution gives the same answers within ~2 %):

| Test | Result |
|---|---|
| Flat ground | 0 iterations, the open wind comes back exactly in every cell |
| 400 m Gaussian hill, wind from the west | Crest +19 % (23.8 km/h), upwind and lee feet -8 %, flanks turned out round the hill by the same amount each side, mass error 19 x smaller |
| 600 m deep north-south valley, wind from the SW | On the floor the wind slows to 16 km/h and turns from 45 deg to 29 deg, towards the valley axis |
| Stability | Flank turning: stable > neutral > unstable |
| Fire model | A field equal to the open wind gives identical spread rates; a field blowing north makes the fire run north |

The mass error is measured from the third layer of elements up. In the lowest two layers the node smoothing
(step 5) mixes in the ground nodes' slip velocity, so div(u) there reflects the smoothing rather than the solve.
The solve itself conserves mass in the finite-element sense everywhere.

### Limits

- **No momentum.** A mass-conserving solver has no momentum equation, so it cannot make flow separation,
  lee-side eddies or reversed flow. In the lee the wind comes out mirroring the windward side, too fast and
  too steady. In field comparisons the WindNinja mass solver did as well as a CFD solver at ridgetops but worse
  on lee slopes (Wagenbrenner et al. 2019, Atmosphere 10: 672). Lee-side fire spread and spotting will be off.
- **Neutral and mechanical only.** There are no daytime upslope or night-time downslope (thermal) winds, which
  dominate when the open wind is light. WindNinja's diurnal option would be the next step.
- **Solved once.** The wind is worked out once for the fixed weather. Changing the wind mid-run would mean
  solving again and recomputing every cell's spread rates while fire arrival times are already queued.
- **Benchmark.** The frame benchmark (`wf_bench [frames]`) still runs on the single open wind. The solver has
  its own benchmark and checksum, `wf_bench --wind` (see Build), and its own tests.

### References

- Forthofer, Butler, Wagenbrenner (2014). A comparison of three approaches for simulating fine-scale surface
  winds in support of wildland fire management. Part I. Int. J. Wildland Fire 23: 969-981.
- Sherman (1978). A mass-consistent model for wind fields over complex terrain. J. Appl. Meteorol. 17: 312-319.
- Wagenbrenner, Forthofer, Page, Butler (2019). Development and evaluation of a Reynolds-averaged Navier-Stokes
  solver in WindNinja for operational wildland fire applications. Atmosphere 10: 672 (in the project files).
- Albini & Baughman (1979), via Andrews (2012) RMRS-GTR-266: the z0 = 0.13 h and d = 0.64 h ratios.

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
| wind_solver | Terrain-adjusted 20-ft wind per cell: mass-conserving finite-element solve on a 100 m terrain-following mesh (WindNinja's method) |
| wind_overlay | Which wind arrows to draw for the current view: map-fixed lattice, length by speed, colour against the open wind (`W`) |
| distance_transform | exact Euclidean distance transform; optionally range-limited to a window round the fire so the per-frame cost follows the fire size. The live fire's field is incremental (only cells that caught fire since the last frame are added) and reads "very far" beyond `LIVE_EDT_RANGE`, the furthest any reader looks |
| safe_planner | fire-safe legs: line of sight, A*, string pulling; reports "no safe path" instead of guessing |
| perimeter | committed fire + sectional perimeter updates; traces the perimeter line into ordered loops |
| voronoi, boustrophedon, aco | the three planning algorithms (boustrophedon: sweep lanes at right angles to the fire, gap fill, ownership map) |
| mission_planner | zones, Voronoi regions (ACO points of interest), sweep lanes and sectors, fleet sizing, sector-to-drone matching, ACO tour |
| drone, mission_controller | drone movement and safety; re-planning, re-routing blocked drones, route assignment; adds and retires drones |
| renderer | all SFML drawing |
| main.cpp | window, events, frame loop |

Library layering: wf_core (everything except GDAL/SFML) <- wf_io (GDAL), wf_render (SFML) <- wildfire.

## Key settings (config.h)

The grid is GRID_WIDTH x GRID_HEIGHT = 1200 x 1000 cells of SIM_CELL_SIZE_M = 10 m (12 km x 10 km). The drone and perimeter
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
of interest per Voronoi region. Only the outer zone is split into Voronoi regions; the inner band is not
segmented (the sweep lanes and sectors cover it).

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

**ACO tour.** The tour only changes when the perimeter does. While the fire grows inside the perimeter (no
section has moved), a re-plan keeps the same Voronoi regions, points of interest and tour order, and only rebuilds
the fire-safe legs against the live fire, so the route stays exactly the same unless the fire has grown across a
leg. The sweep side still re-plans as usual. In the frame benchmark (600 simulated minutes) the perimeter moved 13
times and the tour changed on 14 of 120 re-plans (it used to change on every one).

The tour is always flown anticlockwise. After a re-plan each ACO drone locks back on at the
nearest node of the new tour (or the next one if it has just passed it) and carries on the same way, so drones
never turn back. New drones join in the middle of the largest gaps; when the fleet shrinks, the drone in the most
crowded spot goes home, keeping the rest as evenly spaced as possible. Because drones carry on from where they
were, a new tour can leave them bunched up, with several passing the same node close together. They spread back
out by slowing down: a drone closer than an even share of the tour (tour length / ACO count) behind the drone ahead
flies slower, in proportion to the gap and never below a quarter of full speed (`ACO_SPACING_MIN_SPEED`), until the
gaps are even again; then all fly at full speed. In the tests this takes a few simulated minutes.

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
