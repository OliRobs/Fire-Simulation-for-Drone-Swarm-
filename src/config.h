#pragma once

// --- Simulation grid ---
// The fire runs on GRID_WIDTH x GRID_HEIGHT cells of SIM_CELL_SIZE_M metres: 2000 x 10 m = 20 km square, the
// area of four USGS 1 m DEM tiles (2 x 2 tiles of 10 km). Terrain (elevation, slope, aspect) is worked out
// from the 1 m DEM and averaged into each 10 m cell; the 30 m LANDFIRE fuel and canopy layers are resampled
// onto the same cells. 1 m cells would be 400 million cells (~16 GB) and are finer than the Rothermel model
// is meant for (a uniform fuel bed), so the DEM detail goes into the slope instead.
// The window is centred on the DEM (or the LANDFIRE raster without one) unless --centre <lon> <lat> is given.
constexpr int GRID_WIDTH = 1200;
constexpr int GRID_HEIGHT = 1000;
constexpr int CELL_SIZE = 1;            // Pixel size per grid cell
constexpr float SIM_CELL_SIZE_M = 10.0f;
constexpr float DEFAULT_CELL_SIZE_M = SIM_CELL_SIZE_M; // Ground size of a cell when no DEM is loaded (metres)

// Converts a ground distance into whole cells. The drone and perimeter settings below are design values in
// METRES (they were chosen on 30 m cells), so they stay the same on the ground whatever SIM_CELL_SIZE_M is.
constexpr int cellsFromMetres(float metres) { return static_cast<int>(metres / SIM_CELL_SIZE_M + 0.5f); }

constexpr int PERIMETER_DISTANCE = cellsFromMetres(1200.0f);  // Perimeter ring sits 1.2 km out from the fire
// Lay the window out north-up (local UTM zone) rather than on the LANDFIRE Albers grid, which is rotated away
// from north outside the central US (13 deg in Southern California). Also makes WIND_FROM_DEG a true bearing.
constexpr bool ALIGN_NORTH_UP = true;

// --- Fire spread (Rothermel surface fire model, see rothermel.h) ---
// Simulated time per frame (10 s). With DRONE_SPEED_MPS below, drones fly 16.5 m/s, a realistic survey-drone
// speed, so fire and drones move at the right speeds relative to each other.
constexpr float SIM_MINUTES_PER_TICK = 20.0f / 60.0f;
// Fuel model for cells with no fuel map (synthetic terrain, tests): GR2, low load dry climate grass
constexpr int DEFAULT_FUEL_MODEL = 102;
// Weather. NOT in the LANDFIRE data: these should come from a weather station or forecast (README,
// "Data the Rothermel model still needs").
constexpr float WIND_SPEED_20FT_KMH = 20.0f; // Open wind 20 ft (6.1 m) above the vegetation
constexpr float WIND_FROM_DEG = 225.0f;      // Direction the wind blows FROM, degrees clockwise from north (SW)
constexpr float WIND_ADJUSTMENT_FACTOR = 0.0f; // Midflame / 20-ft wind. <= 0: unsheltered WAF from fuel depth
// Cap the effective wind at Rothermel's "maximum reliable wind" (RMRS-GTR-371 section 3.2.7). On by default:
// with it on, all 40 fuel models reproduce the spread-rate classes Scott & Burgan (2005) publish for them.
constexpr bool APPLY_WIND_LIMIT = true;
// Crown fire from the LANDFIRE canopy layers (crown_fire.h). Foliar moisture is NOT in LANDFIRE: 100 % is the
// usual default (Scott & Reinhardt 2001); RMRS-GTR-220 gives 75-150 % as the normal range.
constexpr bool ENABLE_CROWN_FIRE = true;
constexpr float FOLIAR_MOISTURE_PCT = 100.0f;

// --- LANDFIRE Operational Roads (LF2023_Roads) ---
// Raster values: 20 primary, 21 secondary, 22 tertiary, 23 thinned road; 0 background (from the layer's
// .vat.dbf). Road cells are treated as fuel breaks: no fuel, so surface fire can't spread across them.
constexpr bool ROADS_ARE_FUEL_BREAKS = true;
constexpr int ROAD_CLASS_MIN = 20;
constexpr int ROAD_CLASS_MAX = 23;

// --- Drone survey configuration ---
// Distances are set in metres and converted to cells (cellsFromMetres), so they keep their meaning on the
// ground if the cell size changes. The values are the original 30 m-cell settings (e.g. 70 cells x 30 m).
constexpr int SEARCH_DISTANCE = cellsFromMetres(2100.0f);       // Search area extends 2.1 km from the fire
constexpr int MIN_FIRE_DISTANCE = cellsFromMetres(300.0f);      // Survey standoff: sweeps and points of interest stay 300 m from the fire
constexpr int HARD_MIN_FIRE_DISTANCE = cellsFromMetres(90.0f);  // Absolute limit: flight legs may cut inside the standoff down to 90 m
constexpr double NO_FLY_SQ = (MIN_FIRE_DISTANCE + 0.5) * (MIN_FIRE_DISTANCE + 0.5);           // squared, half-cell margin
constexpr double HARD_NO_FLY_SQ = (HARD_MIN_FIRE_DISTANCE + 0.5) * (HARD_MIN_FIRE_DISTANCE + 0.5);
constexpr int PERIMETER_UPDATE_DISTANCE = cellsFromMetres(600.0f); // A perimeter section moves only once the fire is 600 m from it
constexpr int PERIMETER_SECTION_RADIUS = cellsFromMetres(300.0f);  // How far along the front one update reaches (300 m)
// Range of the live-fire distance transform, in cells. Every reader of the live distance compares it with
// MIN_FIRE_DISTANCE (+0.5) or the smaller HARD_MIN_FIRE_DISTANCE, so distances further out than this are never
// needed and the transform only works out a window this far round the fire (distance_transform.h).
// Derived from the settings above, not a design value of its own.
constexpr int LIVE_EDT_RANGE = MIN_FIRE_DISTANCE + 1;
constexpr int SWEEP_DRONES = 3;           // Boustrophedon drones at the start (the fleet then follows the fire)
constexpr int ACO_DRONES = 5;             // ACO-routed drones at the start (perimeter out to SEARCH_DISTANCE)

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

// --- Shrinking the fleet (both roles) ---
// Growing takes effect at once; shrinking waits. A drone is removed only when one fewer drone would still keep
// every sweep lap (built and timed) or the ACO revisit time within SHRINK_MARGIN of its limit, on SHRINK_REPLANS
// re-plans in a row (re-plans are at least REPLAN_INTERVAL = 10 simulated minutes apart, so 30 minutes or more).
// The removed drone flies a fire-safe route back to the launch point and is taken out of the fleet there.
constexpr float SHRINK_MARGIN = 0.85f;
constexpr int SHRINK_REPLANS = 3;
constexpr float SENSOR_RADIUS = 60.0f / SIM_CELL_SIZE_M;   // Cells; 60 m, so the swath is 120 m wide
constexpr int LANE_SPACING = static_cast<int>(2.0f * SENSOR_RADIUS);
constexpr float DRONE_SPEED_MPS = 16.5f;  // Drone ground speed, m/s
constexpr float DRONE_SPEED = DRONE_SPEED_MPS * SIM_MINUTES_PER_TICK * 60.0f / SIM_CELL_SIZE_M;   // Cells per tick
constexpr int REPLAN_INTERVAL = 60;       // Minimum frames between re-plans as the fire grows
// Voronoi region size target: 630 m x 630 m (was 450 cells of 30 m, ~21 x 21 cells)
constexpr int TARGET_REGION_CELLS = static_cast<int>(405000.0f / (SIM_CELL_SIZE_M * SIM_CELL_SIZE_M) + 0.5f);
constexpr int MIN_SITES = 4;
constexpr int MAX_SITES = 64;
constexpr int LLOYD_ITERATIONS = 4;       // Relaxation passes so regions are evenly sized

// Ant colony optimisation (alpha is fixed at 1, beta at 3)
constexpr int ACO_ANTS = 20;
constexpr int ACO_ITERATIONS = 40;
constexpr float ACO_EVAPORATION = 0.5f;
constexpr float ACO_ELITE_WEIGHT = 2.0f;  // Extra pheromone on the best-so-far tour

// 8-neighbourhood, indexed 0..7
constexpr int NUM_DIRS = 8;
constexpr int DX[NUM_DIRS] = { -1, 0, 1, -1, 1, -1, 0, 1 };
constexpr int DY[NUM_DIRS] = { -1, -1, -1, 0, 0, 1, 1, 1 };
