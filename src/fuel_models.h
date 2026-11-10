#pragma once

// ============================================================================
// The 40 Scott & Burgan (2005) standard fire behaviour fuel models (FBFM40)
//
// LANDFIRE's FBFM40 raster holds a fuel-model NUMBER per 30 m pixel (91 = NB1 urban, 101 = GR1 grass, ...).
// Those numbers are categories, not quantities, so they are looked up here and never averaged or clamped.
// Each entry is one row of Table 7 of Scott & Burgan (2005), RMRS-GTR-153 (in the project docs): the fuel
// inputs Rothermel's surface spread model needs.
// ============================================================================

struct FuelModel {
    int code;              // Fuel model number used in the FBFM40 raster, e.g. 101
    const char* name;      // Fuel model code, e.g. "GR1"
    // Oven-dry fuel load by class, tons/acre (Table 7)
    float load1h, load10h, load100h, loadLiveHerb, loadLiveWoody;
    bool dynamic;          // Dynamic models move live herb load to dead as the herbs cure (Burgan 1979)
    // Surface-area-to-volume ratio, 1/ft (Table 7). 10-h and 100-h are fixed for the whole set: 109 and 30.
    float sav1h, savLiveHerb, savLiveWoody;
    float depthFt;         // Fuel bed depth, ft
    float mxDeadPct;       // Dead fuel moisture of extinction, percent
    float heatBtuLb;       // Heat content of live and dead fuel, Btu/lb
};

// Properties that Scott & Burgan hold constant for every model in the set (RMRS-GTR-153 p. 2 and p. 17)
constexpr float SAV_10H = 109.0f;             // 1/ft
constexpr float SAV_100H = 30.0f;             // 1/ft
constexpr float TOTAL_MINERAL = 0.0555f;      // S_T, fraction
constexpr float EFFECTIVE_MINERAL = 0.010f;   // S_e (silica-free), fraction
constexpr float PARTICLE_DENSITY = 32.0f;     // rho_p, lb/ft^3

// True for the nonburnable block NB1-NB9 (91-99) and for NoData / "no fuel" (<= 0)
bool isNonburnableCode(int code);

// The Table 7 entry for an FBFM40 code, or nullptr if the code is nonburnable or not one of the 40 standard
// models (unknown codes print a warning once, so a custom fuel model in the raster isn't silently ignored).
const FuelModel* findFuelModel(int code);

// All 40 burnable models, for tests and tools
const FuelModel* allFuelModels(int& count);

// Original fuel model 10 (timber litter and understory, Anderson 1982). Not in FBFM40 maps: Rothermel's (1991)
// crown fire spread rate is defined as 3.34 x the surface spread rate of this model (crown_fire.cpp).
const FuelModel& fuelModel10();
