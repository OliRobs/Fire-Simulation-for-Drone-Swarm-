#include "rothermel.h"
#include <algorithm>
#include <cmath>
#include <vector>

// ============================================================================
// Rothermel Surface Fire Spread Model
//
// The basic equation (RMRS-GTR-371 eq. in section 3.1, Table 6):
//
//          I_R * xi * (1 + phi_w + phi_s)        heat source (Btu/ft^2/min)
//     R = --------------------------------  =  -------------------------------
//               rho_b * epsilon * Q_ig           heat sink   (Btu/ft^3)
//
// I_R: reaction intensity, xi: propagating flux ratio, phi_w / phi_s: wind and slope factors,
// rho_b: bulk density, epsilon: effective heating number, Q_ig: heat of preignition.
//
// The fuel bed is a mix of dead (1-h, 10-h, 100-h, cured herb) and live (herb, woody) particle classes,
// combined with the surface-area weighting of Rothermel (1972) and Albini (1976). Comments give the
// symbol names used in RMRS-GTR-371 so each line can be checked against Table 6.
// ============================================================================

namespace {
    // One particle class (one column of Table 7) inside the dead or live category
    struct FuelClass {
        double load;     // w_0, oven-dry load, lb/ft^2
        double sav;      // sigma_ij, 1/ft
        double moisture; // M_f, fraction
    };

    // Albini's (1976) size subclasses for the net-load weighting g_ij (RMRS-GTR-371 Table 6a)
    int sizeSubclass(double sav) {
        if (sav >= 1200.0) return 0;
        if (sav >= 192.0) return 1;
        if (sav >= 96.0) return 2;
        if (sav >= 48.0) return 3;
        if (sav >= 16.0) return 4;
        return 5;        // > 3 inch diameter: ignored (g_ij = 0)
    }

    // Values for one category (dead or live) after surface-area weighting
    struct Category {
        double area = 0.0;          // A_i, mean total surface area per unit fuel cell
        double netLoad = 0.0;       // (w_n)_i = sum_j g_ij (w_n)_ij
        double moisture = 0.0;      // (M_f)_i = sum_j f_ij (M_f)_ij
        double sav = 0.0;           // sigma_i = sum_j f_ij sigma_ij
        double heatSinkTerm = 0.0;  // sum_j f_ij exp(-138/sigma_ij) Q_ig,ij
    };

    // Applies the f_ij and g_ij weightings to one category's classes (only classes with load are passed)
    Category weighCategory(const std::vector<FuelClass>& classes) {
        Category c;
        std::vector<double> area(classes.size());
        for (size_t j = 0; j < classes.size(); ++j) {
            area[j] = classes[j].sav * classes[j].load / PARTICLE_DENSITY;   // A_ij = sigma w_0 / rho_p
            c.area += area[j];
        }
        if (c.area <= 0.0) return c;

        // g_ij: the summed f_ij of every class in the same size subclass
        double subclassF[6] = { 0, 0, 0, 0, 0, 0 };
        for (size_t j = 0; j < classes.size(); ++j) subclassF[sizeSubclass(classes[j].sav)] += area[j] / c.area;

        for (size_t j = 0; j < classes.size(); ++j) {
            const FuelClass& k = classes[j];
            const double f = area[j] / c.area;                                // f_ij = A_ij / A_i
            const int sub = sizeSubclass(k.sav);
            const double g = sub == 5 ? 0.0 : subclassF[sub];                 // g_ij
            c.netLoad += g * k.load * (1.0 - TOTAL_MINERAL);                  // (w_n)_ij = w_0 (1 - S_T)
            c.moisture += f * k.moisture;
            c.sav += f * k.sav;
            const double qig = 250.0 + 1116.0 * k.moisture;                   // Q_ig = 250 + 1116 M_f
            c.heatSinkTerm += f * std::exp(-138.0 / k.sav) * qig;             // epsilon = exp(-138/sigma)
        }
        return c;
    }

    // Moisture damping coefficient eta_M = 1 - 2.59 r + 5.11 r^2 - 3.52 r^3, r = M_f / M_x (max 1)
    double moistureDamping(double mf, double mx) {
        if (mx <= 0.0 || mf >= mx) return 0.0;   // At or above extinction: no combustion (the polynomial is
                                                 // exactly 0 at r = 1, but rounding can leave a tiny residue)
        const double r = mf / mx;
        return std::max(1.0 - 2.59 * r + 5.11 * r * r - 3.52 * r * r * r, 0.0);
    }
}

// ----------------------------------------------------------------------------
// No-wind, no-slope spread for one fuel model
// ----------------------------------------------------------------------------
SurfaceFuelBed computeFuelBed(const FuelModel& m, const FuelMoisture& mo) {
    // --- Dynamic fuel model: move cured herbaceous load from live to dead ---
    // Fraction transferred T = -1.11 M_herb + 1.33, clamped to 0..1 (Burgan 1979; RMRS-GTR-371 Table 7).
    // That is all of it at 30 % herb moisture or less, none at 120 % or more (RMRS-GTR-153 p. 6).
    // The new dead herb class keeps the herb SAV and takes the dead 1-h moisture.
    const double transfer = m.dynamic ? std::clamp(-1.11 * mo.liveHerb + 1.33, 0.0, 1.0) : 0.0;
    const double deadHerb = m.loadLiveHerb * transfer;
    const double liveHerb = m.loadLiveHerb - deadHerb;

    // --- Particle classes (only those carrying load take part) ---
    std::vector<FuelClass> dead, live;
    auto add = [](std::vector<FuelClass>& v, double tonsPerAcre, double sav, double moisture) {
        if (tonsPerAcre > 0.0) v.push_back({ tonsPerAcre * TONS_PER_ACRE_TO_LB_PER_FT2, sav, moisture });
    };
    add(dead, m.load1h, m.sav1h, mo.dead1h);
    add(dead, m.load10h, SAV_10H, mo.dead10h);
    add(dead, m.load100h, SAV_100H, mo.dead100h);
    add(dead, deadHerb, m.savLiveHerb, mo.dead1h);
    add(live, liveHerb, m.savLiveHerb, mo.liveHerb);
    add(live, m.loadLiveWoody, m.savLiveWoody, mo.liveWoody);

    SurfaceFuelBed bed;
    if (dead.empty() || m.depthFt <= 0.0f) return bed;   // Live fuel only burns if dead fuel is present

    const Category d = weighCategory(dead);
    const Category l = weighCategory(live);
    const double totalArea = d.area + l.area;                       // A_T
    const double fDead = d.area / totalArea, fLive = l.area / totalArea;   // f_i = A_i / A_T

    // --- Moisture of extinction ---
    const double mxDead = m.mxDeadPct / 100.0;                     // (M_x)_1, a fuel model input
    double mxLive = mxDead;                                         // (M_x)_2, only matters with live load
    if (!live.empty()) {
        // W: dead-to-live load ratio, M_f,dead: "fine" dead moisture (both only used for (M_x)_2)
        double deadFine = 0.0, deadFineMoist = 0.0, liveFine = 0.0;
        for (const FuelClass& k : dead) {
            const double w = k.load * std::exp(-138.0 / k.sav);
            deadFine += w;
            deadFineMoist += w * k.moisture;
        }
        for (const FuelClass& k : live) liveFine += k.load * std::exp(-500.0 / k.sav);
        const double W = liveFine > 0.0 ? deadFine / liveFine : 0.0;
        const double mfDeadFine = deadFine > 0.0 ? deadFineMoist / deadFine : 0.0;
        mxLive = std::max(2.9 * W * (1.0 - mfDeadFine / mxDead) - 0.226, mxDead);
    }

    // --- Fuel bed properties ---
    double totalLoad = 0.0;
    for (const FuelClass& k : dead) totalLoad += k.load;
    for (const FuelClass& k : live) totalLoad += k.load;
    const double depth = m.depthFt;
    const double rhoB = totalLoad / depth;                          // rho_b, lb/ft^3
    const double beta = rhoB / PARTICLE_DENSITY;                    // beta (one particle density for all)
    const double sigma = fDead * d.sav + fLive * l.sav;             // sigma = sum_i f_i sigma_i
    const double betaOp = 3.348 * std::pow(sigma, -0.8189);         // beta_op
    const double rel = beta / betaOp;

    // --- Reaction intensity I_R = Gamma' sum_i (w_n)_i h_i (eta_M)_i (eta_s)_i ---
    const double s15 = std::pow(sigma, 1.5);
    const double gammaMax = s15 / (495.0 + 0.0594 * s15);          // Gamma'_max
    const double A = 133.0 * std::pow(sigma, -0.7913);              // Albini (1976) form of A
    const double gamma = gammaMax * std::pow(rel, A) * std::exp(A * (1.0 - rel));   // Gamma'
    const double etaS = std::min(0.174 * std::pow(EFFECTIVE_MINERAL, -0.19), 1.0);   // eta_s (max 1)
    const double h = m.heatBtuLb;
    const double iR = gamma * h * etaS * (d.netLoad * moistureDamping(d.moisture, mxDead)
                                        + l.netLoad * moistureDamping(l.moisture, mxLive));

    // --- Propagating flux ratio xi ---
    const double xi = std::exp((0.792 + 0.681 * std::sqrt(sigma)) * (beta + 0.1)) / (192.0 + 0.2595 * sigma);

    // --- Heat sink rho_b sum_i f_i sum_j f_ij epsilon_ij Q_ig,ij ---
    const double heatSink = rhoB * (fDead * d.heatSinkTerm + fLive * l.heatSinkTerm);

    bed.reactionIntensity = std::max(iR, 0.0);
    bed.r0 = heatSink > 0.0 ? bed.reactionIntensity * xi / heatSink : 0.0;
    bed.residenceTime = 384.0 / sigma;                              // t_r, Anderson (1969)
    bed.packingRatio = beta;
    bed.relativePacking = rel;
    bed.sigma = sigma;
    bed.windC = 7.47 * std::exp(-0.133 * std::pow(sigma, 0.55));
    bed.windB = 0.02526 * std::pow(sigma, 0.54);
    bed.windE = 0.715 * std::exp(-3.59e-4 * sigma);
    return bed;
}

// ----------------------------------------------------------------------------
// Wind and slope factors
// ----------------------------------------------------------------------------
double SurfaceFuelBed::windFactor(double u) const {
    if (u <= 0.0) return 0.0;
    return windC * std::pow(u, windB) * std::pow(relativePacking, -windE);
}

double SurfaceFuelBed::slopeFactor(double tanSlope) const {
    if (packingRatio <= 0.0) return 0.0;
    return 5.275 * std::pow(packingRatio, -0.3) * tanSlope * tanSlope;
}

// Inverts the wind factor: U_E = [phi_E (beta/beta_op)^E / C]^(1/B)  (RMRS-GTR-371 section 4.1)
double SurfaceFuelBed::effectiveWindSpeed(double phiE) const {
    if (phiE <= 0.0 || windC <= 0.0 || windB <= 0.0) return 0.0;
    return std::pow(phiE * std::pow(relativePacking, windE) / windC, 1.0 / windB);
}

double unshelteredWaf(double depthFt) {
    if (depthFt <= 0.0) return 0.0;
    return 1.83 / std::log((20.0 + 0.36 * depthFt) / (0.13 * depthFt));
}

double windAdjustmentFactor(double fuelDepthFt, double canopyCoverPct, double canopyHeightFt, double crownRatio) {
    const double f = (canopyCoverPct / 100.0) * std::clamp(crownRatio, 0.0, 1.0) / 3.0;   // Crown fill portion
    if (f <= 0.05 || canopyHeightFt <= 0.0) return unshelteredWaf(fuelDepthFt);         // Open canopy: unsheltered
    const double h = canopyHeightFt;
    return 0.555 / (std::sqrt(f * h) * std::log((20.0 + 0.36 * h) / (0.13 * h)));
}

// ----------------------------------------------------------------------------
// Direction of maximum spread and the fire ellipse
// ----------------------------------------------------------------------------
// The basic equation assumes wind blowing straight upslope. For any other wind direction the wind and
// slope contributions are added as vectors (RMRS-GTR-371 section 6.1, fig. 47):
//     slope vector D_S = R_0 phi_s (upslope), wind vector D_W = R_0 phi_w (downwind)
//     R_H = R_0 + |D_S + D_W|, pointing along D_S + D_W.
// The effective wind factor phi_E = R_H / R_0 - 1 gives the effective wind speed U_E, which sets the
// ellipse's length-to-width ratio Z = 1 + 0.25 U_E (U_E in mi/h) and so its eccentricity.
SpreadEllipse computeSpread(const SurfaceFuelBed& bed, double midflame, double windX, double windY,
                            double tanSlope, double upX, double upY, bool applyWindLimit) {
    SpreadEllipse e;
    if (!bed.burnable()) return e;

    const double phiW = bed.windFactor(midflame);
    const double phiS = bed.slopeFactor(tanSlope);

    // Resultant of the slope and wind vectors (per unit R_0)
    const double vx = phiS * upX + phiW * windX;
    const double vy = phiS * upY + phiW * windY;
    double phiE = std::sqrt(vx * vx + vy * vy);
    if (phiE > 1e-12) { e.dirX = vx / phiE; e.dirY = vy / phiE; }

    double uE = bed.effectiveWindSpeed(phiE);                       // ft/min
    if (applyWindLimit) {
        // Maximum reliable wind U = 0.9 I_R (Rothermel 1972 p. 33); spread stops increasing beyond it
        const double uMax = 0.9 * bed.reactionIntensity;
        if (uE > uMax) {
            uE = uMax;
            phiE = bed.windFactor(uE);
        }
    }

    e.headRos = bed.r0 * (1.0 + phiE);
    e.effectiveWindMph = uE / FT_PER_MIN_PER_MPH;
    const double z = 1.0 + 0.25 * e.effectiveWindMph;               // Length-to-width ratio
    e.eccentricity = std::sqrt(z * z - 1.0) / z;
    return e;
}

double SpreadEllipse::rosToward(double ux, double uy) const {
    const double cosGamma = ux * dirX + uy * dirY;                  // Angle from the head direction
    return headRos * (1.0 - eccentricity) / (1.0 - eccentricity * cosGamma);
}
