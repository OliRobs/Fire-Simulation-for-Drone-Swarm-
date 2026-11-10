#include "fuel_models.h"
#include <iostream>
#include <set>

// ============================================================================
// FBFM40 Fuel-Model Table
//
// Source: Scott, J.H. & Burgan, R.E. (2005) Standard Fire Behavior Fuel Models: A Comprehensive Set for Use
// with Rothermel's Surface Fire Spread Model. USDA Forest Service RMRS-GTR-153, Table 7 (p. 18).
// Values are copied exactly as published (US customary units); rothermel.cpp converts them.
//
// Where Table 7 lists an SAV of 9999 the class has no load, so the value never enters the calculation;
// it is kept here only so the rows read the same as the paper.
// ============================================================================

namespace {
    constexpr bool DYN = true;   // Dynamic: herbaceous load is split between live and dead by curing
    constexpr bool STA = false;  // Static, or no live herbaceous load ("N/A" in Table 7)

    constexpr FuelModel kModels[] = {
        //  code  name    1-h   10-h  100-h  herb  woody  type  SAV1h herb  woody depth Mx  heat
        // --- Grass (GR) ---
        { 101, "GR1", 0.10f, 0.00f, 0.00f, 0.30f, 0.00f, DYN, 2200, 2000, 9999, 0.4f, 15, 8000 },
        { 102, "GR2", 0.10f, 0.00f, 0.00f, 1.00f, 0.00f, DYN, 2000, 1800, 9999, 1.0f, 15, 8000 },
        { 103, "GR3", 0.10f, 0.40f, 0.00f, 1.50f, 0.00f, DYN, 1500, 1300, 9999, 2.0f, 30, 8000 },
        { 104, "GR4", 0.25f, 0.00f, 0.00f, 1.90f, 0.00f, DYN, 2000, 1800, 9999, 2.0f, 15, 8000 },
        { 105, "GR5", 0.40f, 0.00f, 0.00f, 2.50f, 0.00f, DYN, 1800, 1600, 9999, 1.5f, 40, 8000 },
        { 106, "GR6", 0.10f, 0.00f, 0.00f, 3.40f, 0.00f, DYN, 2200, 2000, 9999, 1.5f, 40, 9000 },
        { 107, "GR7", 1.00f, 0.00f, 0.00f, 5.40f, 0.00f, DYN, 2000, 1800, 9999, 3.0f, 15, 8000 },
        { 108, "GR8", 0.50f, 1.00f, 0.00f, 7.30f, 0.00f, DYN, 1500, 1300, 9999, 4.0f, 30, 8000 },
        { 109, "GR9", 1.00f, 1.00f, 0.00f, 9.00f, 0.00f, DYN, 1800, 1600, 9999, 5.0f, 40, 8000 },
        // --- Grass-shrub (GS) ---
        { 121, "GS1", 0.20f, 0.00f, 0.00f, 0.50f, 0.65f, DYN, 2000, 1800, 1800, 0.9f, 15, 8000 },
        { 122, "GS2", 0.50f, 0.50f, 0.00f, 0.60f, 1.00f, DYN, 2000, 1800, 1800, 1.5f, 15, 8000 },
        { 123, "GS3", 0.30f, 0.25f, 0.00f, 1.45f, 1.25f, DYN, 1800, 1600, 1600, 1.8f, 40, 8000 },
        { 124, "GS4", 1.90f, 0.30f, 0.10f, 3.40f, 7.10f, DYN, 1800, 1600, 1600, 2.1f, 40, 8000 },
        // --- Shrub (SH) ---
        { 141, "SH1", 0.25f, 0.25f, 0.00f, 0.15f, 1.30f, DYN, 2000, 1800, 1600, 1.0f, 15, 8000 },
        { 142, "SH2", 1.35f, 2.40f, 0.75f, 0.00f, 3.85f, STA, 2000, 9999, 1600, 1.0f, 15, 8000 },
        { 143, "SH3", 0.45f, 3.00f, 0.00f, 0.00f, 6.20f, STA, 1600, 9999, 1400, 2.4f, 40, 8000 },
        { 144, "SH4", 0.85f, 1.15f, 0.20f, 0.00f, 2.55f, STA, 2000, 1800, 1600, 3.0f, 30, 8000 },
        { 145, "SH5", 3.60f, 2.10f, 0.00f, 0.00f, 2.90f, STA,  750, 9999, 1600, 6.0f, 15, 8000 },
        { 146, "SH6", 2.90f, 1.45f, 0.00f, 0.00f, 1.40f, STA,  750, 9999, 1600, 2.0f, 30, 8000 },
        { 147, "SH7", 3.50f, 5.30f, 2.20f, 0.00f, 3.40f, STA,  750, 9999, 1600, 6.0f, 15, 8000 },
        { 148, "SH8", 2.05f, 3.40f, 0.85f, 0.00f, 4.35f, STA,  750, 9999, 1600, 3.0f, 40, 8000 },
        { 149, "SH9", 4.50f, 2.45f, 0.00f, 1.55f, 7.00f, DYN,  750, 1800, 1500, 4.4f, 40, 8000 },
        // --- Timber-understory (TU) ---
        { 161, "TU1", 0.20f, 0.90f, 1.50f, 0.20f, 0.90f, DYN, 2000, 1800, 1600, 0.6f, 20, 8000 },
        { 162, "TU2", 0.95f, 1.80f, 1.25f, 0.00f, 0.20f, STA, 2000, 9999, 1600, 1.0f, 30, 8000 },
        { 163, "TU3", 1.10f, 0.15f, 0.25f, 0.65f, 1.10f, DYN, 1800, 1600, 1400, 1.3f, 30, 8000 },
        { 164, "TU4", 4.50f, 0.00f, 0.00f, 0.00f, 2.00f, STA, 2300, 9999, 2000, 0.5f, 12, 8000 },
        { 165, "TU5", 4.00f, 4.00f, 3.00f, 0.00f, 3.00f, STA, 1500, 9999,  750, 1.0f, 25, 8000 },
        // --- Timber litter (TL) ---
        { 181, "TL1", 1.00f, 2.20f, 3.60f, 0.00f, 0.00f, STA, 2000, 9999, 9999, 0.2f, 30, 8000 },
        { 182, "TL2", 1.40f, 2.30f, 2.20f, 0.00f, 0.00f, STA, 2000, 9999, 9999, 0.2f, 25, 8000 },
        { 183, "TL3", 0.50f, 2.20f, 2.80f, 0.00f, 0.00f, STA, 2000, 9999, 9999, 0.3f, 20, 8000 },
        { 184, "TL4", 0.50f, 1.50f, 4.20f, 0.00f, 0.00f, STA, 2000, 9999, 9999, 0.4f, 25, 8000 },
        { 185, "TL5", 1.15f, 2.50f, 4.40f, 0.00f, 0.00f, STA, 2000, 9999, 1600, 0.6f, 25, 8000 },
        { 186, "TL6", 2.40f, 1.20f, 1.20f, 0.00f, 0.00f, STA, 2000, 9999, 9999, 0.3f, 25, 8000 },
        { 187, "TL7", 0.30f, 1.40f, 8.10f, 0.00f, 0.00f, STA, 2000, 9999, 9999, 0.4f, 25, 8000 },
        { 188, "TL8", 5.80f, 1.40f, 1.10f, 0.00f, 0.00f, STA, 1800, 9999, 9999, 0.3f, 35, 8000 },
        { 189, "TL9", 6.65f, 3.30f, 4.15f, 0.00f, 0.00f, STA, 1800, 9999, 1600, 0.6f, 35, 8000 },
        // --- Slash-blowdown (SB) ---
        { 201, "SB1", 1.50f, 3.00f, 11.00f, 0.00f, 0.00f, STA, 2000, 9999, 9999, 1.0f, 25, 8000 },
        { 202, "SB2", 4.50f, 4.25f, 4.00f, 0.00f, 0.00f, STA, 2000, 9999, 9999, 1.0f, 25, 8000 },
        { 203, "SB3", 5.50f, 2.75f, 3.00f, 0.00f, 0.00f, STA, 2000, 9999, 9999, 1.2f, 25, 8000 },
        { 204, "SB4", 5.25f, 3.50f, 5.25f, 0.00f, 0.00f, STA, 2000, 9999, 9999, 2.7f, 25, 8000 },
    };
    constexpr int kModelCount = static_cast<int>(sizeof(kModels) / sizeof(kModels[0]));
}

// ----------------------------------------------------------------------------
// Nonburnable codes
// ----------------------------------------------------------------------------
// NB1 urban (91), NB2 snow/ice (92), NB3 agriculture (93), NB8 water (98), NB9 bare ground (99), plus the
// rest of the reserved 90-99 block. Codes <= 0 cover LANDFIRE's -9999 NoData and "no fuel loaded".
bool isNonburnableCode(int code) {
    return code <= 0 || (code >= 90 && code <= 99);
}

// ----------------------------------------------------------------------------
// Code -> Table 7 row
// ----------------------------------------------------------------------------
const FuelModel* findFuelModel(int code) {
    if (isNonburnableCode(code)) return nullptr;

    for (const FuelModel& m : kModels) {
        if (m.code == code) return &m;
    }

    // Custom or unexpected code: warn once per code. It is then treated as nonburnable rather than guessed,
    // because a made-up fuel model would give made-up spread rates.
    static std::set<int> warned;
    if (warned.insert(code).second) {
        std::cerr << "[Fuel] FBFM40 code " << code << " is not one of the 40 standard models; "
                     "treating it as nonburnable\n";
    }
    return nullptr;
}

const FuelModel* allFuelModels(int& count) {
    count = kModelCount;
    return kModels;
}
