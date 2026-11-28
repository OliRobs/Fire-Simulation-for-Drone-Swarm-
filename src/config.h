#pragma once

// ============================================================================
// Settings. In five sections:
//   1. Grid             the app's default grid and the cell size (the code takes a GridSpec, grid.h)
//   2. Fire weather     DEFAULTS only: a scenario file replaces them at start-up (wildfire --weather, weather.h)
//   3. Model settings   fire, terrain-wind and road handling choices
//   4. Design limits    the drone-survey behaviour: distances, fleet sizes, time limits, speeds
//   5. Derived values   worked out from the above; not settings of their own
// plus the neighbour tables at the end. The values in sections 3 and 4 are design decisions: ask the owner
// before changing them (CLAUDE.md).
// ============================================================================

// ============================================================================
// 1. Grid
// ============================================================================
// The fire runs on GRID_WIDTH x GRID_HEIGHT cells of SIM_CELL_SIZE_M metres: 1200 x 1000 cells of 10 m =
// 12 km east-west by 10 km north-south. That is a window cut from the four USGS 1 m DEM tiles (2 x 2 tiles of
// 10 km, 20 km square together), centred on them, so the tiles reach about 4 km past the window east and west
// and 5 km north and south. Terrain (elevation, slope, aspect) is worked out from the 1 m DEM and averaged into
// each 10 m cell; the 30 m LANDFIRE fuel and canopy layers are resampled onto the same cells. 1 m cells would
// be 120 million cells for this window (~5 GB) and are finer than the Rothermel model is meant for (a uniform
// fuel bed), so the DEM detail goes into the slope instead.
// The window is centred on the DEM (or the LANDFIRE raster without one) unless --centre <lon> <lat> is given.
// GRID_WIDTH x GRID_HEIGHT is only the app's default size: main.cpp and wf_bench build a GridSpec from it and
// everything else takes that (grid.h).
constexpr int GRID_WIDTH = 1200;
constexpr int GRID_HEIGHT = 1000;
constexpr int CELL_SIZE = 1;            // Pixel size per grid cell
constexpr float SIM_CELL_SIZE_M = 10.0f;
constexpr float DEFAULT_CELL_SIZE_M = SIM_CELL_SIZE_M; // Ground size of a cell when no DEM is loaded (metres)

// Converts a ground distance into whole cells. The drone and perimeter settings below are design values in
// METRES (they were chosen on 30 m cells), so they stay the same on the ground whatever SIM_CELL_SIZE_M is.
constexpr int cellsFromMetres(float metres) { return static_cast<int>(metres / SIM_CELL_SIZE_M + 0.5f); }

// Lay the window out north-up (local UTM zone) rather than on the LANDFIRE Albers grid, which is rotated away
// from north outside the central US (13 deg in Southern California). Also makes WIND_FROM_DEG a true bearing.
constexpr bool ALIGN_NORTH_UP = true;

// Simulated time per frame (5 s). With DRONE_SPEED_MPS below, drones fly 16.5 m/s, a realistic survey-drone
// speed, so fire and drones move at the right speeds relative to each other.
constexpr float SIM_MINUTES_PER_TICK = 5.0f / 60.0f;

// ============================================================================
// 2. Fire weather: DEFAULTS (FireWeather, defaultWeather() in weather.h)
// ============================================================================
// NOT in the LANDFIRE data: these should come from a weather station or forecast (README, "Data the Rothermel
// model still needs"). A scenario file replaces any of them at start-up without a rebuild
// (wildfire --weather <file>; scenarios/default.ini holds these same values with their sources), and
// --wind / --wind-from override the wind. The fuel moisture defaults are in rothermel.h (FuelMoisture).
constexpr float WIND_SPEED_20FT_KMH = 20.0f; // Open wind 20 ft (6.1 m) above the vegetation
constexpr float WIND_FROM_DEG = 225.0f;      // Direction the wind blows FROM, degrees clockwise from north (SW)
constexpr float WIND_ADJUSTMENT_FACTOR = 0.0f; // Midflame / 20-ft wind. <= 0: unsheltered WAF from fuel depth
// Cap the effective wind at Rothermel's "maximum reliable wind" (RMRS-GTR-371 section 3.2.7). On by default:
// with it on, all 40 fuel models reproduce the spread-rate classes Scott & Burgan (2005) publish for them.
constexpr bool APPLY_WIND_LIMIT = true;
// Atmospheric stability for the terrain wind as alpha = alpha_H / alpha_V: 1 neutral (as WindNinja without its
// stability option), > 1 unstable (more air goes over hills), < 1 stable (more air goes round them).
// WindNinja's table runs from 0.2 (Pasquill class F, clear calm night) to 5 (class A, strong sun, light wind).
constexpr float WIND_STABILITY_ALPHA = 1.0f;
// Crown fire from the LANDFIRE canopy layers (crown_fire.h). Foliar moisture is NOT in LANDFIRE: 100 % is the
// usual default (Scott & Reinhardt 2001); RMRS-GTR-220 gives 75-150 % as the normal range.
constexpr bool ENABLE_CROWN_FIRE = true;
constexpr float FOLIAR_MOISTURE_PCT = 100.0f;

// ============================================================================
// 3. Model settings
// ============================================================================
// --- Fire spread (Rothermel surface fire model, see rothermel.h) ---
// Fuel model for cells with no fuel map (synthetic terrain, tests): GR2, low load dry climate grass
constexpr int DEFAULT_FUEL_MODEL = 102;

// --- Terrain-adjusted wind (wind_solver.h) ---
// The open wind above is adjusted over the terrain by a mass-conserving solver (WindNinja's method), so every
// cell gets its own 20-ft wind speed and direction: faster over ridges, slower in hollows, turned along valleys.
// false = the single open wind everywhere, as before.
constexpr bool USE_TERRAIN_WIND = true;
constexpr float WIND_COLUMN_SIZE_M = 100.0f;   // Horizontal spacing of the wind mesh (10 x 10 fire cells per column)
constexpr int WIND_LAYERS = 20;                // Mesh nodes up each column, ground to domain top (WindNinja: 20)
// Height of the wind mesh's top above the highest summit. The top is an open boundary, so a low top lets the air
// escape upwards instead of speeding up over the ridges. WindNinja's rule (column width / 400, grown over the
// layers) puts it only ~160 m up at 100 m columns, which cut the crest speed-up of a 400 m test hill from ~20 %
// to ~2 %; the result stops changing once the top is ~2.5 km up (wind_solver.h).
constexpr float WIND_TOP_ABOVE_SUMMIT_M = 2500.0f;

// --- LANDFIRE Operational Roads (LF2023_Roads) ---
// Raster values: 20 primary, 21 secondary, 22 tertiary, 23 thinned road; 0 background (from the layer's
// .vat.dbf). Road cells are treated as fuel breaks: no fuel, so surface fire can't spread across them.
constexpr bool ROADS_ARE_FUEL_BREAKS = true;
constexpr int ROAD_CLASS_MIN = 20;
constexpr int ROAD_CLASS_MAX = 23;

// ============================================================================
// 4. Design limits: the drone survey
// ============================================================================
// --- Distances round the fire ---
// Set in metres and converted to cells (cellsFromMetres), so they keep their meaning on the ground if the cell
// size changes. The values are the original 30 m-cell settings (e.g. 70 cells x 30 m).
constexpr int PERIMETER_DISTANCE = cellsFromMetres(1200.0f);  // Perimeter ring sits 1.2 km out from the fire
constexpr int SEARCH_DISTANCE = cellsFromMetres(2100.0f);       // Search area extends 2.1 km from the fire
constexpr int MIN_FIRE_DISTANCE = cellsFromMetres(300.0f);      // Survey standoff: sweeps and points of interest stay 300 m from the fire
constexpr int HARD_MIN_FIRE_DISTANCE = cellsFromMetres(90.0f);  // Absolute limit: flight legs may cut inside the standoff down to 90 m
constexpr int PERIMETER_UPDATE_DISTANCE = cellsFromMetres(600.0f); // A perimeter section moves only once the fire is 600 m from it
constexpr int PERIMETER_SECTION_RADIUS = cellsFromMetres(300.0f);  // How far along the front one update reaches (300 m)

// --- Drones ---
constexpr int SWEEP_DRONES = 3;           // Boustrophedon drones at the start (the fleet then follows the fire)
constexpr int ACO_DRONES = 2;             // ACO-routed drones at the start (perimeter out to SEARCH_DISTANCE)
constexpr float DRONE_SPEED_MPS = 16.5f;  // Drone ground speed, m/s
constexpr int REPLAN_INTERVAL = 60;       // Minimum frames between re-plans as the fire grows

// --- Sweep fleet sizing ---
// Each sweep drone owns one sector of the inner band (a contiguous run of lanes along the perimeter line) and
// flies it as a closed loop. Lap time = length of that closed route (lanes, steps between lanes, fire detours
// and the closing leg back to the start) / DRONE_SPEED_MPS. Because the route is a closed loop, the lap time
// is also the longest wait between visits to any lane. The fleet is sized so no lap is over the limit.
constexpr float MAX_SWEEP_LAP_S = 20.0f * 60.0f;  // Longest allowed sweep lap: 20 minutes
constexpr float SWEEP_LAP_MARGIN = 0.9f;          // Sectors are first sized to 90 % of the limit (straight-line
                                                  // estimate), leaving 10 % for fire detours
constexpr int MIN_SWEEP_DRONES = 3;               // Fewest sweep drones, even for a small fire
constexpr int MAX_SWEEP_DRONES = 10;              // Most sweep drones; above this, sectors may go over the limit

// --- ACO fleet sizing ---
// The ACO drones share one closed tour, evenly spaced, each flying the whole tour. With n drones every point on
// the tour is passed every (one-drone loop time / n), so n = ceil(loop time / ACO_REVISIT_S). Adding drones
// shortens the wait between visits, not any one drone's loop.
constexpr float ACO_REVISIT_S = 20.0f * 60.0f;    // Longest allowed wait between visits to a point on the tour
constexpr int MIN_ACO_DRONES = 1;                 // Fewest ACO drones (one flies the whole tour early on)
constexpr int MAX_ACO_DRONES = 10;                // Most ACO drones; above this, the revisit time may be exceeded
// Keeping ACO drones evenly spaced. Drones never turn back, and can't fly faster than DRONE_SPEED, so a drone
// closer than one even share of the tour (tour length / ACO count) behind the drone ahead flies slower, in
// proportion to the gap, until the gap has opened back up. Within ACO_SPACING_TOLERANCE of an even share it
// flies at full speed, so small differences from waypoint rounding don't keep slowing drones down.
constexpr float ACO_SPACING_MIN_SPEED = 0.25f;    // Slowest a crowded ACO drone flies, as a fraction of DRONE_SPEED
constexpr float ACO_SPACING_TOLERANCE = 0.05f;    // Gaps within 5 % of an even share count as evenly spaced

// --- Shrinking the fleet (both roles) ---
// Growing takes effect at once; shrinking waits. A drone is removed only when one fewer drone would still keep
// every sweep lap (built and timed) or the ACO revisit time within SHRINK_MARGIN of its limit, on SHRINK_REPLANS
// re-plans in a row (re-plans are at least REPLAN_INTERVAL = 10 simulated minutes apart, so 30 minutes or more).
// The removed drone flies a fire-safe route back to the launch point and is taken out of the fleet there.
constexpr float SHRINK_MARGIN = 0.85f;
constexpr int SHRINK_REPLANS = 3;

// --- Voronoi regions (points of interest for the ACO tour) ---
constexpr int MIN_SITES = 4;
constexpr int MAX_SITES = 64;
constexpr int LLOYD_ITERATIONS = 4;       // Relaxation passes so regions are evenly sized

// --- Ant colony optimisation (alpha is fixed at 1, beta at 3) ---
constexpr int ACO_ANTS = 20;
constexpr int ACO_ITERATIONS = 40;
constexpr float ACO_EVAPORATION = 0.5f;
constexpr float ACO_ELITE_WEIGHT = 2.0f;  // Extra pheromone on the best-so-far tour

// ============================================================================
// 5. Derived values (worked out from the settings above)
// ============================================================================
constexpr double NO_FLY_SQ = (MIN_FIRE_DISTANCE + 0.5) * (MIN_FIRE_DISTANCE + 0.5);           // squared, half-cell margin
constexpr double HARD_NO_FLY_SQ = (HARD_MIN_FIRE_DISTANCE + 0.5) * (HARD_MIN_FIRE_DISTANCE + 0.5);
// Range of the live-fire distance transform, in cells. Every reader of the live distance compares it with
// MIN_FIRE_DISTANCE (+0.5) or the smaller HARD_MIN_FIRE_DISTANCE, so distances further out than this are never
// needed and the transform only works out a window this far round the fire (distance_transform.h).
constexpr int LIVE_EDT_RANGE = MIN_FIRE_DISTANCE + 1;
constexpr float SENSOR_RADIUS = 60.0f / SIM_CELL_SIZE_M;   // Cells; 60 m, so the swath is 120 m wide
constexpr int LANE_SPACING = static_cast<int>(2.0f * SENSOR_RADIUS);
constexpr float DRONE_SPEED = DRONE_SPEED_MPS * SIM_MINUTES_PER_TICK * 60.0f / SIM_CELL_SIZE_M;   // Cells per tick
// Voronoi region size target: 630 m x 630 m (was 450 cells of 30 m, ~21 x 21 cells)
constexpr int TARGET_REGION_CELLS = static_cast<int>(405000.0f / (SIM_CELL_SIZE_M * SIM_CELL_SIZE_M) + 0.5f);

// ============================================================================
// 8-neighbourhood, indexed 0..7
// ============================================================================
constexpr int NUM_DIRS = 8;
constexpr int DX[NUM_DIRS] = { -1, 0, 1, -1, 1, -1, 0, 1 };
constexpr int DY[NUM_DIRS] = { -1, -1, -1, 0, 0, 1, 1, 1 };
