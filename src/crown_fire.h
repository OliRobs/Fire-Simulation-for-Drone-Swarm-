#pragma once
#include <cstdint>
#include "rothermel.h"

// ============================================================================
// Crown fire: when a surface fire climbs into the forest canopy, and how fast it then spreads
//
// Uses the LANDFIRE canopy layers (canopy base height CBH, canopy bulk density CBD) with the method FARSITE
// and FlamMap use by default (Finney 1998, "Finney" crown fire method in RMRS-GTR-220, the LANDFIRE guidebook
// in the project docs):
//   1. Initiation (Van Wagner 1977): the surface fire torches the canopy once its fireline intensity reaches
//          I_0 = (0.010 * CBH * (460 + 25.9 * FMC))^1.5    kW/m, CBH in m, FMC = foliar moisture in %
//   2. Crown spread rate (Rothermel 1991): R_crown = 3.34 x the surface spread rate of original fuel model 10,
//      with the midflame wind taken as 0.4 x the 20-ft wind.
//   3. Active or passive (Van Wagner 1977): crown fire carries itself from crown to crown only if
//          R_crown >= R'_active = 3.0 / CBD    m/min, CBD in kg/m^3
//      Otherwise it is passive (single trees or groups torch) and spreads at the surface rate.
// ============================================================================

enum class FireType : std::uint8_t {
    Surface = 0,   // Fire stays on the ground
    Passive = 1,   // Torching: trees catch, but the fire still moves at the surface rate
    Active = 2     // Running crown fire at the crown spread rate
};

// Van Wagner (1977) critical surface fireline intensity for crown ignition, kW/m
double criticalSurfaceIntensity(double canopyBaseHeightM, double foliarMoisturePct);

// Van Wagner (1977) critical (minimum) spread rate for active crown fire, m/min
double criticalActiveRos(double canopyBulkDensityKgM3);

// Byram's fireline intensity I_B = H_A R, H_A = I_R t_r (RMRS-GTR-371 sections 4.3-4.4), converted to kW/m.
// rosFtPerMin: the spread rate in the direction of interest (normally the head fire).
double firelineIntensityKwM(const SurfaceFuelBed& bed, double rosFtPerMin);

// Byram's flame length L = 0.0775 I_B^0.46 (m, I_B in kW/m), for reporting and the guidebook's rule of thumb
double flameLengthM(double firelineIntensityKwM);
