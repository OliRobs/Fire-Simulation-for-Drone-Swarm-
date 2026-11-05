#pragma once

// --- Configuration & Constants ---
constexpr int GRID_WIDTH = 500;
constexpr int GRID_HEIGHT = 500;
constexpr int CELL_SIZE = 2;            // Pixel size per grid cell
constexpr int PERIMETER_DISTANCE = 40;  // Perimeter ring sits this many cells out from the fire
constexpr float BASE_IGNITION_RATE = 0.005f;
constexpr float FUEL_BURN_RATE = 0.01f;

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
