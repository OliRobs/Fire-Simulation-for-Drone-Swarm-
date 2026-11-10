#pragma once
#include "fuel_models.h"

// ============================================================================
// Rothermel (1972) surface fire spread model, with the Albini (1976) revisions
//
// Equations as listed in Andrews, P.L. (2018) The Rothermel surface fire spread model and associated
// developments: a comprehensive explanation. USDA Forest Service RMRS-GTR-371 (in the project docs):
// Table 6 (final model for heterogeneous fuel), section 4 (residence time, effective wind speed) and
// section 6 / Table 26 (direction of maximum spread and the elliptical fire shape, as in BehavePlus).
//
// Everything here works in the paper's US customary units (ft, lb, Btu, min). Callers convert at the edges:
// spread rates come back in ft/min, and FT_PER_MIN_TO_M_PER_MIN converts them.
// No SFML, no GDAL: this is part of wf_core.
// ============================================================================

constexpr double FT_PER_M = 3.280839895;
constexpr double FT_PER_MIN_TO_M_PER_MIN = 0.3048;
constexpr double KMH_TO_FT_PER_MIN = 1000.0 * FT_PER_M / 60.0;  // 1 km/h = 54.68 ft/min
constexpr double FT_PER_MIN_PER_MPH = 88.0;                     // 1 mi/h = 88 ft/min
constexpr double TONS_PER_ACRE_TO_LB_PER_FT2 = 2000.0 / 43560.0; // 0.0459137 (RMRS-GTR-371 Table 27 note)

// Fuel moisture contents as FRACTIONS of oven-dry weight (0.06 = 6 %).
// These are NOT in LANDFIRE: they come from weather (see README "Data the Rothermel model still needs").
// The defaults are Scott & Burgan's standard scenario D2L2 (RMRS-GTR-153 Tables 3 and 4): low dead fuel
// moisture, two-thirds cured herbaceous, low live woody moisture.
struct FuelMoisture {
    double dead1h = 0.06;
    double dead10h = 0.07;
    double dead100h = 0.08;
    double liveHerb = 0.60;
    double liveWoody = 0.90;
};

// No-wind, no-slope result for one fuel model and moisture, plus the fuel-bed values that the wind and
// slope factors need. Computed once per fuel model (it doesn't depend on the cell).
struct SurfaceFuelBed {
    double r0 = 0.0;                // No-wind no-slope rate of spread, ft/min (0 = will not spread)
    double reactionIntensity = 0.0; // I_R, Btu/ft^2/min
    double residenceTime = 0.0;     // t_r = 384 / sigma, min (Anderson 1969)
    double packingRatio = 0.0;      // beta
    double relativePacking = 0.0;   // beta / beta_op
    double sigma = 0.0;             // Characteristic surface-area-to-volume ratio, 1/ft
    double windC = 0.0, windB = 0.0, windE = 0.0;  // Wind factor coefficients C, B, E

    bool burnable() const { return r0 > 0.0; }

    // Wind factor phi_w = C U^B (beta/beta_op)^-E, U = midflame wind speed in ft/min
    double windFactor(double midflameFtPerMin) const;
    // Slope factor phi_s = 5.275 beta^-0.3 tan^2(phi)
    double slopeFactor(double tanSlope) const;
    // Effective midflame wind speed (ft/min) that would give the effective wind factor phiE on its own
    double effectiveWindSpeed(double phiE) const;
};

// Rothermel no-wind no-slope calculation for a Scott & Burgan fuel model, including the dynamic transfer of
// cured herbaceous load from live to dead (Burgan 1979) for dynamic models.
SurfaceFuelBed computeFuelBed(const FuelModel& model, const FuelMoisture& moisture);

// Unsheltered wind adjustment factor (midflame / 20-ft wind) from fuel bed depth in ft:
// WAF = 1.83 / ln((20 + 0.36 H) / (0.13 H)) (Albini & Baughman 1979, as given in Andrews 2012, RMRS-GTR-266,
// which RMRS-GTR-371 refers to). "Unsheltered" means no forest canopy above the surface fuel.
double unshelteredWaf(double fuelDepthFt);

// Wind adjustment factor allowing for a forest canopy (Andrews 2012, RMRS-GTR-266, after Albini & Baughman 1979).
// Crown fill portion f = (canopy cover / 100) x crown ratio / 3 (crowns taken as cones). If f <= 0.05 the fuel is
// unsheltered and unshelteredWaf(depth) is used; otherwise the fuel is sheltered by the canopy:
//     WAF = 0.555 / ( sqrt(f H) ln((20 + 0.36 H) / (0.13 H)) ),  H = canopy height in ft.
// crownRatio = (canopy height - canopy base height) / canopy height, 0..1.
double windAdjustmentFactor(double fuelDepthFt, double canopyCoverPct, double canopyHeightFt, double crownRatio);

// Head fire rate and direction for one cell, from vector addition of the wind and slope effects
// (RMRS-GTR-371 section 6.1, Table 26), and the ellipse that gives the rate in any other direction.
// Directions are unit vectors in any consistent 2D frame (the simulation uses grid x = east, y = south).
struct SpreadEllipse {
    double headRos = 0.0;            // R_H, ft/min, in the direction of maximum spread
    double dirX = 1.0, dirY = 0.0;   // Unit vector of maximum spread (arbitrary when there's no wind/slope)
    double eccentricity = 0.0;       // e = sqrt(Z^2 - 1) / Z, Z = 1 + 0.25 U_E (U_E in mi/h)
    double effectiveWindMph = 0.0;   // U_E, combined effect of wind and slope

    // Rate of spread from the ignition point in direction (ux, uy) (unit vector):
    // R_gamma = R_H (1 - e) / (1 - e cos gamma)
    double rosToward(double ux, double uy) const;
    double backingRos() const { return headRos * (1.0 - eccentricity) / (1.0 + eccentricity); }
};

// windDirX/Y: unit vector the wind blows TOWARDS. upslopeX/Y: unit vector pointing uphill (ignored when
// tanSlope is 0). applyWindLimit: cap the effective wind at the "maximum reliable wind" U = 0.9 I_R
// (Rothermel 1972, as in BehavePlus); RMRS-GTR-371 section 3.2.7 discusses why this cap is often left off.
SpreadEllipse computeSpread(const SurfaceFuelBed& bed, double midflameFtPerMin, double windDirX, double windDirY,
                            double tanSlope, double upslopeX, double upslopeY, bool applyWindLimit);
