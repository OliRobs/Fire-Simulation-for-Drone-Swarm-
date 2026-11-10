#include "fuel_models.h"
#include <iostream>
#include <set>

// ============================================================================
// FBFM40 Fuel-Model Lookup
//
// Source: Scott, J.H. & Burgan, R.E. (2005) Standard Fire Behavior Fuel Models,
// USDA Forest Service RMRS-GTR-153 (in the project docs).
//
//  - Spread class: the adjective spread-rate class each model is given in the fuel model selection guide
//    (pp. 10-12), for 5 mi/h midflame wind, dry dead fuel, two-thirds cured herbaceous, zero slope.
//  - Total load: the sum of 1-h + 10-h + 100-h + live herb + live woody load from Table 7 (tons/acre).
//
// Mapping onto the automaton (these two scalings are our own calibration choice, not from the paper):
//  - flammability = spread class mapped to 0.25 (very low) ... 1.0 (extreme)
//  - load         = 0.2 + 0.8 * totalLoad / max total load, so heavy shrub and slash burn far longer than grass
// ============================================================================

namespace {
    // Spread-rate adjective classes from the selection guide
    enum Spread { VL, L, M, H, VH, E };

    constexpr float kFlammability[] = { 0.25f, 0.40f, 0.60f, 0.80f, 0.90f, 1.00f };

    struct Entry {
        int code;           // FBFM40 fuel model number
        Spread spread;      // Spread-rate class (selection guide)
        float totalLoad;    // Total fuel load, tons/acre (Table 7)
    };

    // All 40 burnable standard fuel models
    constexpr Entry kModels[] = {
        // Grass (GR1-GR9)
        { 101, M, 0.40f }, { 102, H, 1.10f }, { 103, H, 2.00f }, { 104, VH, 2.15f }, { 105, VH, 2.90f },
        { 106, VH, 3.50f }, { 107, VH, 6.40f }, { 108, VH, 8.80f }, { 109, E, 11.00f },
        // Grass-shrub (GS1-GS4)
        { 121, M, 1.35f }, { 122, H, 2.60f }, { 123, H, 3.25f }, { 124, H, 12.80f },
        // Shrub (SH1-SH9)
        { 141, VL, 1.95f }, { 142, L, 8.35f }, { 143, L, 9.65f }, { 144, H, 4.75f }, { 145, VH, 8.60f },
        { 146, H, 5.75f }, { 147, H, 14.40f }, { 148, H, 10.65f }, { 149, H, 15.50f },
        // Timber-understory (TU1-TU5)
        { 161, L, 3.70f }, { 162, M, 4.20f }, { 163, H, 3.25f }, { 164, M, 6.50f }, { 165, M, 14.00f },
        // Timber litter (TL1-TL9)
        { 181, VL, 6.80f }, { 182, VL, 5.90f }, { 183, VL, 5.50f }, { 184, L, 6.20f }, { 185, L, 8.05f },
        { 186, M, 4.80f }, { 187, L, 9.80f }, { 188, M, 8.30f }, { 189, M, 14.10f },
        // Slash-blowdown (SB1-SB4)
        { 201, M, 15.50f }, { 202, M, 12.75f }, { 203, H, 11.25f }, { 204, VH, 14.00f },
    };

    constexpr float kMaxTotalLoad = 15.50f;  // SH9 and SB1, the heaviest models in the set

    FuelProperties fromEntry(const Entry& e) {
        return { kFlammability[e.spread], 0.2f + 0.8f * e.totalLoad / kMaxTotalLoad };
    }
}

FuelProperties fuelFromFbfm40(int code) {
    // NoData (LANDFIRE uses -9999, filled to 0 by the loader) and the non-burnable block NB1-NB9:
    // urban, snow/ice, agriculture, open water, bare ground. No fuel, so fire cannot spread there.
    if (code <= 0 || (code >= 90 && code <= 99)) return { 0.0f, 0.0f };

    for (const Entry& e : kModels) {
        if (e.code == code) return fromEntry(e);
    }

    // Custom or unexpected code: warn once per code, then treat it as a moderate fuel
    static std::set<int> warned;
    if (warned.insert(code).second) {
        std::cerr << "[Fuel] Unknown FBFM40 code " << code << ", treating as moderate fuel\n";
    }
    return { kFlammability[M], 0.5f };
}
