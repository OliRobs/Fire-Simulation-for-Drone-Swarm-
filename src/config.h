#pragma once

// --- Configuration & Constants ---
constexpr int GRID_WIDTH = 500;
constexpr int GRID_HEIGHT = 500;
constexpr int CELL_SIZE = 2;            // Pixel size per grid cell
constexpr int PERIMETER_DISTANCE = 40;  // Perimeter ring sits this many cells out from the fire
constexpr float DEFAULT_CELL_SIZE_M = 30.0f; // Ground size of a cell when no DEM is loaded (metres)

// --- Simulation window on the LANDFIRE rasters ---
// LANDFIRE layers are 30 m pixels covering a large area, so the simulation runs on a GRID_WIDTH x GRID_HEIGHT
// window of them (500 x 30 m = 15 km square). The window is centred on the raster's
// centre unless --centre <lon> <lat> is given on the command line.
constexpr float SIM_CELL_SIZE_M = 30.0f;
// Lay the window out north-up (local UTM zone) rather than on the LANDFIRE Albers grid, which is rotated away
// from north outside the central US (13 deg in Southern California). Also makes WIND_FROM_DEG a true bearing.
constexpr bool ALIGN_NORTH_UP = true;

// --- Fire spread (Rothermel surface fire model, see rothermel.h) ---
// Simulated time per frame. With 30 m cells, DRONE_SPEED (5.5 cells per tick) is then 16.5 m/s, a realistic
// survey-drone speed, so fire and drones move at the right speeds relative to each other.
constexpr float SIM_MINUTES_PER_TICK = 10.0f / 60.0f;
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
constexpr int SEARCH_DISTANCE = 70;       // Search area extends this many cells from the fire
constexpr int MIN_FIRE_DISTANCE = 10;      // Survey standoff: sweeps and points of interest stay this far from the fire
constexpr int HARD_MIN_FIRE_DISTANCE = 3; // Absolute limit: flight legs may cut inside the standoff down to this
constexpr double NO_FLY_SQ = (MIN_FIRE_DISTANCE + 0.5) * (MIN_FIRE_DISTANCE + 0.5);           // squared, half-cell margin
constexpr double HARD_NO_FLY_SQ = (HARD_MIN_FIRE_DISTANCE + 0.5) * (HARD_MIN_FIRE_DISTANCE + 0.5);
constexpr int PERIMETER_UPDATE_DISTANCE = 20;  // A perimeter section moves only once the fire is this close to it
constexpr int PERIMETER_SECTION_RADIUS = 10;  // How far along the front one update reaches (cells)
constexpr int SWEEP_DRONES = 3;           // Boustrophedon drones (inside the perimeter)
constexpr int ACO_DRONES = 5;             // ACO-routed drones (perimeter out to SEARCH_DISTANCE)
constexpr float SENSOR_RADIUS = 2.0f;     // Cells; swath width = 2 * radius
constexpr int LANE_SPACING = static_cast<int>(2.0f * SENSOR_RADIUS);
constexpr float DRONE_SPEED = 5.5f;       // Cells per tick
constexpr int REPLAN_INTERVAL = 60;       // Minimum frames between re-plans as the fire grows
constexpr int TARGET_REGION_CELLS = 450;  // Voronoi region size target (~21 x 21 cells)
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
