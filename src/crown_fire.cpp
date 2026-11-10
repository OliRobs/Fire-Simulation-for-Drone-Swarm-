#include "crown_fire.h"
#include <cmath>
#include <limits>

// ============================================================================
// Crown Fire Thresholds (see crown_fire.h for the method and sources)
// ============================================================================

namespace {
    constexpr double BTU_PER_FT_S_TO_KW_PER_M = 3.46141;   // 1055.06 J / 0.3048 m
}

// ----------------------------------------------------------------------------
// Van Wagner (1977) eq. 4, with his empirical constant C = 0.010
// ----------------------------------------------------------------------------
// The heat needed to raise foliage to ignition is 460 + 25.9 FMC kJ/kg; a higher canopy base or wetter
// foliage needs a hotter surface fire. Example (RMRS-GTR-220 Table 6): CBH 5 m, FMC 100 % gives about
// 1900 kW/m, i.e. a 2.5 m flame.
double criticalSurfaceIntensity(double cbhM, double fmcPct) {
    if (cbhM <= 0.0) return 0.0;
    return std::pow(0.010 * cbhM * (460.0 + 25.9 * fmcPct), 1.5);
}

// ----------------------------------------------------------------------------
// Van Wagner (1977): active crown fire needs a mass flow rate of canopy fuel of at least 3.0 kg/m^2/min,
// so the spread rate must be at least 3.0 / CBD.
// ----------------------------------------------------------------------------
double criticalActiveRos(double cbdKgM3) {
    if (cbdKgM3 <= 0.0) return std::numeric_limits<double>::infinity();   // No canopy fuel: never active
    return 3.0 / cbdKgM3;
}

// ----------------------------------------------------------------------------
// Byram's fireline intensity: heat per unit area (Btu/ft^2) x spread rate (ft/min) / 60 = Btu/ft/s
// ----------------------------------------------------------------------------
double firelineIntensityKwM(const SurfaceFuelBed& bed, double rosFtPerMin) {
    const double heatPerArea = bed.reactionIntensity * bed.residenceTime;   // H_A = I_R t_r
    return heatPerArea * rosFtPerMin / 60.0 * BTU_PER_FT_S_TO_KW_PER_M;
}

// Byram (1959) flame length in SI form
double flameLengthM(double iKwM) {
    return iKwM > 0.0 ? 0.0775 * std::pow(iKwM, 0.46) : 0.0;
}
