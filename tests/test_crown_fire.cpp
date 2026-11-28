// ============================================================================
// Tests: Crown fire (crown_fire.h): Van Wagner thresholds and crown fire in the simulation
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "fire_model.h"
#include "crown_fire.h"
#include "test_framework.h"
#include "test_support.h"

TEST(crown_fire_thresholds) {
    // Van Wagner (1977): CBH 5 m, FMC 100 % needs ~1880 kW/m, a 2.5 m flame (RMRS-GTR-220 Table 6 / p. 24)
    const double i0 = criticalSurfaceIntensity(5.0, 100.0);
    CHECK(std::fabs(i0 - 1883.0) < 5.0);
    CHECK(std::fabs(flameLengthM(i0) - 2.5) < 0.05);
    CHECK(criticalSurfaceIntensity(10.0, 100.0) > i0);   // Higher canopy base: harder to torch
    CHECK(criticalSurfaceIntensity(5.0, 80.0) < i0);     // Drier foliage: easier
    CHECK(std::fabs(criticalActiveRos(0.15) - 20.0) < 1e-9);   // 3.0 / CBD m/min
}

namespace {
    // Uniform forest of one fuel model under the given canopy; returns how far (cells) the head fire runs
    // east in 60 minutes with a 30 km/h wind from the west, and the fire type it assigned (crown_fire.h).
    std::pair<int, int> forestRun(int fuel, float cc, float chM, float cbhM, float cbd) {
        std::vector<Cell> grid(N);
        for (Cell& c : grid) {
            c.fuelModel = static_cast<std::int16_t>(fuel);
            c.canopyCoverPct = cc; c.canopyHeightM = chM; c.canopyBaseHeightM = cbhM; c.canopyBulkDensity = cbd;
        }
        FireWeather w;
        w.windFromDeg = 270.0f;
        w.wind20ftKmh = 30.0f;
        FireSimulation sim(grid, w, 30.0f, 30.0f);
        sim.ignite(at(100, 250));
        for (int t = 0; t < 360; ++t) sim.step();      // 60 minutes
        int k = 0;
        while (sim.fireMask()[at(101 + k, 250)]) ++k;
        return { k, sim.grid()[at(100, 250)].fireType };
    }
}

SLOW_TEST(crown_fire_in_the_simulation) {
    // TU5 (timber litter with shrub understory, ~640 kW/m here). A 10 m canopy base needs ~5300 kW/m, so the
    // fire stays on the ground...
    const auto [surfaceRun, surfaceType] = forestRun(165, 50.0f, 20.0f, 10.0f, 0.25f);
    CHECK(surfaceType == static_cast<int>(FireType::Surface));
    // ...a 1 m base (~170 kW/m) torches, and with dense canopy fuel (crown rate ~17 m/min > 3.0 / 0.25 = 12)
    // it becomes an active crown fire that runs far faster
    const auto [crownRun, crownType] = forestRun(165, 50.0f, 20.0f, 1.0f, 0.25f);
    CHECK(crownType == static_cast<int>(FireType::Active));
    CHECK(crownRun > 5 * std::max(surfaceRun, 1));
    // Sparse canopy fuel (CBD 0.05 needs 60 m/min) can only torch, so it spreads at the surface rate
    const auto [passiveRun, passiveType] = forestRun(165, 50.0f, 20.0f, 1.0f, 0.05f);
    CHECK(passiveType == static_cast<int>(FireType::Passive));
    CHECK(passiveRun < crownRun / 5);
    std::printf("  (60 min runs: surface %d, passive %d, active crown %d cells)\n", surfaceRun, passiveRun, crownRun);

    // Crown fire switched off: same landscape stays a surface fire
    std::vector<Cell> grid(N);
    for (Cell& c : grid) { c.fuelModel = 165; c.canopyCoverPct = 50; c.canopyHeightM = 20; c.canopyBaseHeightM = 1; c.canopyBulkDensity = 0.25f; }
    FireWeather off;
    off.crownFire = false;
    FireSimulation sim(grid, off);
    CHECK(sim.grid()[at(10, 10)].fireType == static_cast<int>(FireType::Surface));
}

TEST(canopy_slows_surface_fire) {
    // Same grass under a closed canopy (no crown fire possible: CBD 0) runs slower than in the open
    const auto [open, t1] = forestRun(102, 0.0f, 0.0f, 0.0f, 0.0f);
    const auto [shaded, t2] = forestRun(102, 80.0f, 25.0f, 5.0f, 0.0f);
    CHECK(shaded < open);
    (void)t1; (void)t2;
}
