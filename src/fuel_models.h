#pragma once

// Converts LANDFIRE FBFM40 fuel-model codes (Scott & Burgan 2005, RMRS-GTR-153) into the two numbers the
// cellular automaton uses. FBFM40 rasters hold category codes (91 = urban, 98 = water, 101 = GR1 grass, ...),
// not densities, so they must be looked up, never clamped or averaged.
struct FuelProperties {
    float flammability; // 0..1, scales the chance that fire spreads INTO the cell (0 = nonburnable)
    float load;         // 0..1, sets how long the cell burns once alight (burn ticks = load / FUEL_BURN_RATE)
};

// Unknown codes print a warning once and are treated as a moderate fuel.
// Non-burnable codes (91-99) and NoData (<= 0) return { 0, 0 }.
FuelProperties fuelFromFbfm40(int code);
