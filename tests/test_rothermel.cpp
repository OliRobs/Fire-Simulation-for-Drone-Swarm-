// ============================================================================
// Tests: Rothermel surface fire model (rothermel.h) and fire-model set-up rates
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "fire_model.h"
#include "crown_fire.h"
#include "fuel_models.h"
#include "rothermel.h"
#include "test_framework.h"
#include "test_support.h"

// ============================================================================
// Fire model
// ============================================================================
// --- Rothermel model ---
namespace {
    // Scott & Burgan's adjective spread-rate classes (RMRS-GTR-153 Table 5), in chains per hour
    enum SpreadClass { VL, L, M, H, VH, E };
    SpreadClass classify(double chainsPerHour) {
        if (chainsPerHour < 2) return VL;
        if (chainsPerHour < 5) return L;
        if (chainsPerHour < 20) return M;
        if (chainsPerHour < 50) return H;
        if (chainsPerHour < 150) return VH;
        return E;
    }
    // Head fire rate (ft/min) for a fuel model on flat ground with a midflame wind blowing along +x
    double headRosFtMin(const FuelModel& m, const FuelMoisture& mo, double midflameMph, bool windLimit) {
        const SurfaceFuelBed bed = computeFuelBed(m, mo);
        return computeSpread(bed, midflameMph * FT_PER_MIN_PER_MPH, 1, 0, 0, 0, 0, windLimit).headRos;
    }
}

// Every fuel model's head fire rate must fall in the spread class the selection guide gives it
// (RMRS-GTR-153 pp. 10-12: moisture scenario D2L2, 5 mi/h midflame wind, zero slope).
TEST(rothermel_reproduces_scott_burgan_spread_classes) {
    const std::pair<const char*, SpreadClass> expected[] = {
        { "GR1", M }, { "GR2", H }, { "GR3", H }, { "GR4", VH }, { "GR5", VH }, { "GR6", VH }, { "GR7", VH },
        { "GR8", VH }, { "GR9", E }, { "GS1", M }, { "GS2", H }, { "GS3", H }, { "GS4", H },
        { "SH1", VL }, { "SH2", L }, { "SH3", L }, { "SH4", H }, { "SH5", VH }, { "SH6", H }, { "SH7", H },
        { "SH8", H }, { "SH9", H }, { "TU1", L }, { "TU2", M }, { "TU3", H }, { "TU4", M }, { "TU5", M },
        { "TL1", VL }, { "TL2", VL }, { "TL3", VL }, { "TL4", L }, { "TL5", L }, { "TL6", M }, { "TL7", L },
        { "TL8", M }, { "TL9", M }, { "SB1", M }, { "SB2", M }, { "SB3", H }, { "SB4", VH },
    };
    int count = 0;
    const FuelModel* models = allFuelModels(count);
    CHECK(count == 40);
    for (int i = 0; i < count; ++i) {
        for (const auto& [name, cls] : expected) {
            if (std::string(name) != models[i].name) continue;
            const double chainsPerHour = headRosFtMin(models[i], FuelMoisture(), 5.0, true) * 60.0 / 66.0;
            if (classify(chainsPerHour) != cls) {
                std::printf("  %s: %.1f ch/h\n", name, chainsPerHour);
                CHECK(classify(chainsPerHour) == cls);
            }
        }
    }
}

// Packing ratio and characteristic SAV against the summary values on the fuel model pages (RMRS-GTR-153)
TEST(rothermel_fuel_bed_matches_published_values) {
    struct Ref { int code; double beta; double sigma; };
    for (const Ref r : { Ref{ 101, 0.00143, 2054 }, Ref{ 142, 0.01198, 1672 }, Ref{ 181, 0.04878, 1716 },
                         Ref{ 164, 0.01865, 2216 }, Ref{ 146, 0.00412, 1144 } }) {
        const SurfaceFuelBed bed = computeFuelBed(*findFuelModel(r.code), FuelMoisture());
        CHECK(std::fabs(bed.packingRatio - r.beta) / r.beta < 0.005);
        CHECK(std::fabs(bed.sigma - r.sigma) < 1.0);
    }
}

TEST(rothermel_moisture_and_curing) {
    const FuelModel& gr2 = *findFuelModel(102);
    // Dead fuel wetter than the 15 % moisture of extinction: no spread at all
    FuelMoisture wet;
    wet.dead1h = wet.dead10h = wet.dead100h = 0.30;
    CHECK(!computeFuelBed(gr2, wet).burnable());
    // Dynamic model: fully cured grass (herb moisture 30 %) burns much faster than green grass (120 %)
    FuelMoisture cured, green;
    cured.liveHerb = 0.30;
    green.liveHerb = 1.20;
    CHECK(computeFuelBed(gr2, cured).r0 > 2.0 * computeFuelBed(gr2, green).r0);
}

TEST(rothermel_wind_slope_and_ellipse) {
    const SurfaceFuelBed bed = computeFuelBed(*findFuelModel(145), FuelMoisture());   // SH5
    // No wind, no slope: a circle, spreading at R_0 every way
    const SpreadEllipse still = computeSpread(bed, 0, 1, 0, 0, 0, 0, false);
    CHECK(std::fabs(still.headRos - bed.r0) < 1e-9 && still.eccentricity == 0.0);
    CHECK(std::fabs(still.rosToward(0, 1) - bed.r0) < 1e-9);

    // Wind towards +x: head fire faster, flanks slower, backing slowest
    const SpreadEllipse windy = computeSpread(bed, 5 * FT_PER_MIN_PER_MPH, 1, 0, 0, 0, 0, false);
    CHECK(windy.headRos > 3.0 * bed.r0);
    CHECK(std::fabs(windy.rosToward(1, 0) - windy.headRos) < 1e-9);
    CHECK(windy.rosToward(0, 1) < windy.headRos && windy.backingRos() < windy.rosToward(0, 1));

    // Slope alone (30 degrees, uphill = -y): head points uphill and gets the slope factor
    const double tan30 = std::tan(30.0 * 3.14159265358979 / 180.0);
    const SpreadEllipse slope = computeSpread(bed, 0, 1, 0, tan30, 0, -1, false);
    CHECK(std::fabs(slope.dirY + 1.0) < 1e-9);
    CHECK(std::fabs(slope.headRos - bed.r0 * (1 + bed.slopeFactor(tan30))) < 1e-9);

    // Wind across the slope: the head direction lies between the two
    const SpreadEllipse both = computeSpread(bed, 5 * FT_PER_MIN_PER_MPH, 1, 0, tan30, 0, -1, false);
    CHECK(both.dirX > 0.0 && both.dirY < 0.0);

    // Unsheltered WAF: deeper fuel beds feel more of the 20-ft wind (Andrews 2012)
    CHECK(unshelteredWaf(1.0) > 0.3 && unshelteredWaf(1.0) < 0.4);
    CHECK(unshelteredWaf(6.0) > unshelteredWaf(1.0));
}

// --- Rothermel set-up: cached fuel-bed terms and the per-cell memos give bit-identical rates ---
namespace {
    // The wind and slope factors and the effective wind speed exactly as they were written before the fuel-bed
    // terms were cached (kept word for word, so the new code can be checked against them bit for bit)
    double origWindFactor(const SurfaceFuelBed& b, double u) {
        if (u <= 0.0) return 0.0;
        return b.windC * std::pow(u, b.windB) * std::pow(b.relativePacking, -b.windE);
    }
    double origSlopeFactor(const SurfaceFuelBed& b, double tanSlope) {
        if (b.packingRatio <= 0.0) return 0.0;
        return 5.275 * std::pow(b.packingRatio, -0.3) * tanSlope * tanSlope;
    }
    double origEffectiveWindSpeed(const SurfaceFuelBed& b, double phiE) {
        if (phiE <= 0.0 || b.windC <= 0.0 || b.windB <= 0.0) return 0.0;
        return std::pow(phiE * std::pow(b.relativePacking, b.windE) / b.windC, 1.0 / b.windB);
    }

    // The original computeSpread, word for word, on the original factors above
    SpreadEllipse origComputeSpread(const SurfaceFuelBed& bed, double midflame, double windX, double windY,
                                    double tanSlope, double upX, double upY, bool applyWindLimit) {
        SpreadEllipse e;
        if (!bed.burnable()) return e;
        const double phiW = origWindFactor(bed, midflame);
        const double phiS = origSlopeFactor(bed, tanSlope);
        const double vx = phiS * upX + phiW * windX;
        const double vy = phiS * upY + phiW * windY;
        double phiE = std::sqrt(vx * vx + vy * vy);
        if (phiE > 1e-12) { e.dirX = vx / phiE; e.dirY = vy / phiE; }
        double uE = origEffectiveWindSpeed(bed, phiE);
        if (applyWindLimit) {
            const double uMax = 0.9 * bed.reactionIntensity;
            if (uE > uMax) {
                uE = uMax;
                phiE = origWindFactor(bed, uE);
            }
        }
        e.headRos = bed.r0 * (1.0 + phiE);
        e.effectiveWindMph = uE / FT_PER_MIN_PER_MPH;
        const double z = 1.0 + 0.25 * e.effectiveWindMph;
        e.eccentricity = std::sqrt(z * z - 1.0) / z;
        return e;
    }

    // Moisture scenarios used by both tests: the D2L2 default, very dry, fully cured, and wet enough that some
    // models stop burning (so non-burnable beds, with all-zero fields, are exercised too)
    std::vector<FuelMoisture> testMoistures() {
        FuelMoisture dry, cured, wet;
        dry.dead1h = 0.03; dry.dead10h = 0.04; dry.dead100h = 0.05; dry.liveHerb = 0.30; dry.liveWoody = 0.60;
        cured.liveHerb = 0.30;
        wet.dead1h = 0.20; wet.dead10h = 0.22; wet.dead100h = 0.24; wet.liveHerb = 1.50; wet.liveWoody = 1.50;
        return { FuelMoisture(), dry, cured, wet };
    }
}

// The wind factor, slope factor, effective wind speed and whole spread ellipse, now using the fuel-bed terms
// cached by computeFuelBed, must equal the original inline-pow formulas bit for bit: every fuel model (plus
// fuel model 10), four moisture scenarios, and a spread of wind speeds, slopes and directions, with and
// without the wind limit. computeSpreadWithWindFactor with a precomputed phi_w must equal computeSpread.
TEST(rothermel_cached_bed_terms_match_original_formulas) {
    int count = 0;
    const FuelModel* models = allFuelModels(count);
    std::vector<const FuelModel*> list;
    for (int i = 0; i < count; ++i) list.push_back(&models[i]);
    list.push_back(&fuelModel10());

    std::mt19937 rng(11);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    long compared = 0, differ = 0;
    int burnableBeds = 0, deadBeds = 0;
    for (const FuelMoisture& mo : testMoistures()) {
        for (const FuelModel* m : list) {
            const SurfaceFuelBed bed = computeFuelBed(*m, mo);
            (bed.burnable() ? burnableBeds : deadBeds)++;
            for (int k = 0; k < 200; ++k) {
                // Wind 0..60 mi/h midflame (k = 0: exactly 0), slope 0..45 degrees, random wind and aspect
                const double u = k == 0 ? 0.0 : 60.0 * FT_PER_MIN_PER_MPH * unit(rng);
                const double tanSlope = k % 7 == 0 ? 0.0 : std::tan(0.785 * unit(rng));
                const double phiE = 50.0 * unit(rng);
                const double wa = 6.283185307 * unit(rng), sa = 6.283185307 * unit(rng);
                const bool limit = k % 2 == 1;
                ++compared;
                if (bed.windFactor(u) != origWindFactor(bed, u)) ++differ;
                if (bed.slopeFactor(tanSlope) != origSlopeFactor(bed, tanSlope)) ++differ;
                if (bed.effectiveWindSpeed(phiE) != origEffectiveWindSpeed(bed, phiE)) ++differ;
                const SpreadEllipse a = computeSpread(bed, u, std::cos(wa), std::sin(wa), tanSlope,
                                                      std::cos(sa), std::sin(sa), limit);
                const SpreadEllipse b = origComputeSpread(bed, u, std::cos(wa), std::sin(wa), tanSlope,
                                                          std::cos(sa), std::sin(sa), limit);
                const SpreadEllipse c = computeSpreadWithWindFactor(bed, bed.windFactor(u), std::cos(wa),
                                                                    std::sin(wa), tanSlope, std::cos(sa),
                                                                    std::sin(sa), limit);
                for (const SpreadEllipse* s : { &a, &c }) {
                    if (s->headRos != b.headRos || s->dirX != b.dirX || s->dirY != b.dirY
                        || s->eccentricity != b.eccentricity || s->effectiveWindMph != b.effectiveWindMph) ++differ;
                }
            }
        }
    }
    std::printf("  %ld cases (%d burnable beds, %d that won't burn), %ld differ\n", compared, burnableBeds,
                deadBeds, differ);
    CHECK(differ == 0);
    CHECK(burnableBeds > 0 && deadBeds > 0);   // Both kinds of bed were exercised
}

// FireSimulation's set-up (fuel-bed terms cached, wind factor memo per fuel model, unsheltered WAF per fuel
// model, crown wind factor computed once) must give exactly the per-cell rates of the original set-up loop,
// kept word for word below on the original formulas. The grid mixes every fuel model in patches, random
// slopes and aspects, forest patches with crown fire (active and passive), awkward canopy edge cases (cover
// without height, height without cover, sparse cover on the 5 % crown fill boundary, canopy with no bulk
// density), roads, non-burnable
// codes and code 0. Several weathers: the default, a fixed WAF, the wind limit, no crown fire, calm, and a
// wet scenario. Cells are 10 x 12 m so the step directions are not the square-grid ones.
SLOW_TEST(fire_setup_matches_original_rates) {
    int count = 0;
    const FuelModel* models = allFuelModels(count);
    std::mt19937 rng(23);
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    std::vector<Cell> grid(N);
    for (int y = 0; y < GRID_HEIGHT; ++y) {
        for (int x = 0; x < GRID_WIDTH; ++x) {
            Cell& c = grid[at(x, y)];
            const int patch = (x / 30) * 7 + (y / 30) * 13;
            c.fuelModel = static_cast<std::int16_t>(models[patch % count].code);
            if (patch % 53 == 0) c.fuelModel = 91;   // Non-burnable (urban)
            if (patch % 61 == 0) c.fuelModel = 0;    // No data
            if (x % 97 == 0) c.roadClass = static_cast<std::uint8_t>(20 + y % 4);   // A road: fuel break
            c.slopeDeg = (x / 30 + y / 30) % 5 == 0 ? 0.0f : 45.0f * unit(rng);   // Some flat patches
            const float a = 6.2831853f * unit(rng);
            c.upslopeX = c.slopeDeg > 0.0f ? std::cos(a) : 0.0f;
            c.upslopeY = c.slopeDeg > 0.0f ? std::sin(a) : 0.0f;
            switch (((x / 40) + 3 * (y / 40)) % 6) {
            case 0: case 1:   // Forest
                c.canopyCoverPct = 10.0f + 80.0f * unit(rng);
                c.canopyHeightM = 5.0f + 30.0f * unit(rng);
                c.canopyBaseHeightM = 0.3f + 6.0f * unit(rng);
                c.canopyBulkDensity = 0.02f + 0.35f * unit(rng);
                break;
            case 2:           // Odd canopy layers: cover without height, height without cover, sparse cover,
                              // and a canopy with no bulk density (shelters the wind, can't crown)
                if (x % 4 == 0) c.canopyCoverPct = 50.0f;
                else if (x % 4 == 1) { c.canopyHeightM = 20.0f; c.canopyBulkDensity = 0.1f; }
                else if (x % 4 == 2) { c.canopyCoverPct = 15.0f * unit(rng); c.canopyHeightM = 15.0f;
                                       c.canopyBaseHeightM = 15.0f * unit(rng); c.canopyBulkDensity = 0.1f; }
                else { c.canopyCoverPct = 60.0f; c.canopyHeightM = 20.0f; c.canopyBaseHeightM = 2.0f; }
                break;
            default: break;   // Open ground: no canopy
            }
        }
    }

    std::vector<FireWeather> weathers;
    weathers.push_back(FireWeather());
    { FireWeather w; w.windAdjustmentFactor = 0.4f; weathers.push_back(w); }
    { FireWeather w; w.windLimit = !w.windLimit; w.wind20ftKmh = 60.0f; weathers.push_back(w); }
    { FireWeather w; w.crownFire = false; w.windFromDeg = 135.0f; weathers.push_back(w); }
    { FireWeather w; w.wind20ftKmh = 0.0f; weathers.push_back(w); }
    { FireWeather w; w.moisture = testMoistures()[3]; w.wind20ftKmh = 45.0f; weathers.push_back(w); }

    const float csx = 10.0f, csy = 12.0f;
    long cellsCompared = 0, differ = 0, active = 0, passive = 0, sheltered = 0;
    for (const FireWeather& weather : weathers) {
        FireSimulation sim(grid, weather, csx, csy);

        // --- The original set-up loop, word for word except for writing to local arrays ---
        double stepUx[NUM_DIRS], stepUy[NUM_DIRS];
        for (int d = 0; d < NUM_DIRS; ++d) {
            const double sx = DX[d] * csx, sy = DY[d] * csy;
            const double len = std::sqrt(sx * sx + sy * sy);
            stepUx[d] = sx / len;
            stepUy[d] = sy / len;
        }
        const Vec2 wind = windTowardGrid(weather.windFromDeg);
        const double wind20 = weather.wind20ftKmh * KMH_TO_FT_PER_MIN;
        const SurfaceFuelBed fm10 = computeFuelBed(fuelModel10(), weather.moisture);
        constexpr double DEG_TO_RAD = 0.017453292519943295;
        const double minCell = std::min(csx, csy);
        for (int idx = 0; idx < N; ++idx) {
            const Cell& c = grid[idx];
            std::uint8_t fireType = static_cast<std::uint8_t>(FireType::Surface);
            bool burnable = false;
            float head = 0.0f, dur = 0.0f, dir[NUM_DIRS] = {};
            const int code = c.fuelModel;
            const FuelModel* m = code > 0 && code < 256 ? findFuelModel(code) : nullptr;
            const SurfaceFuelBed bed = m ? computeFuelBed(*m, weather.moisture) : SurfaceFuelBed();
            const bool road = ROADS_ARE_FUEL_BREAKS && c.roadClass >= ROAD_CLASS_MIN && c.roadClass <= ROAD_CLASS_MAX;
            if (m && bed.burnable() && !road) {
                const double tanSlope = std::tan(c.slopeDeg * DEG_TO_RAD);
                double waf = weather.windAdjustmentFactor;
                if (waf <= 0.0) {
                    const double crownRatio = c.canopyHeightM > 0.0f
                        ? (c.canopyHeightM - c.canopyBaseHeightM) / c.canopyHeightM : 0.0;
                    waf = windAdjustmentFactor(m->depthFt, c.canopyCoverPct, c.canopyHeightM * FT_PER_M, crownRatio);
                    if (waf != unshelteredWaf(m->depthFt)) ++sheltered;
                }
                SpreadEllipse e = origComputeSpread(bed, wind20 * waf, wind.x, wind.y, tanSlope, c.upslopeX,
                                                    c.upslopeY, weather.windLimit);
                if (weather.crownFire && c.canopyBulkDensity > 0.0f && c.canopyHeightM > 0.0f && c.canopyCoverPct > 0.0f) {
                    const double surfaceIntensity = firelineIntensityKwM(bed, e.headRos);
                    if (surfaceIntensity >= criticalSurfaceIntensity(c.canopyBaseHeightM, weather.foliarMoisturePct)) {
                        SpreadEllipse crown = origComputeSpread(fm10, wind20 * 0.4, wind.x, wind.y, tanSlope,
                                                                c.upslopeX, c.upslopeY, false);
                        crown.headRos *= 3.34;
                        const double crownMPerMin = crown.headRos * FT_PER_MIN_TO_M_PER_MIN;
                        if (crownMPerMin >= criticalActiveRos(c.canopyBulkDensity)) {
                            fireType = static_cast<std::uint8_t>(FireType::Active);
                            if (crown.headRos > e.headRos) e = crown;
                        }
                        else {
                            fireType = static_cast<std::uint8_t>(FireType::Passive);
                        }
                    }
                }
                burnable = true;
                head = static_cast<float>(e.headRos * FT_PER_MIN_TO_M_PER_MIN);
                for (int d = 0; d < NUM_DIRS; ++d)
                    dir[d] = static_cast<float>(e.rosToward(stepUx[d], stepUy[d]) * FT_PER_MIN_TO_M_PER_MIN);
                dur = static_cast<float>(minCell / head + bed.residenceTime);
            }

            // --- Compare bit for bit with what FireSimulation set up ---
            ++cellsCompared;
            if (fireType == static_cast<std::uint8_t>(FireType::Active)) ++active;
            if (fireType == static_cast<std::uint8_t>(FireType::Passive)) ++passive;
            bool same = sim.canBurn(idx) == burnable && sim.grid()[idx].fireType == fireType
                        && sim.headRos(idx) == head;
            if (burnable) {
                same = same && sim.burnDuration(idx) == dur;
                for (int d = 0; d < NUM_DIRS; ++d) same = same && sim.directionalRos(idx, d) == dir[d];
            }
            if (!same) ++differ;
        }
    }
    std::printf("  %ld cells over %zu weathers (%ld active crown, %ld passive, %ld sheltered), %ld differ\n",
                cellsCompared, weathers.size(), active, passive, sheltered, differ);
    CHECK(differ == 0);
    CHECK(active > 0 && passive > 0 && sheltered > 0);   // The crown fire and canopy paths were exercised
}

// --- Canopy and crown fire ---
TEST(fuel_model_10_matches_published_values) {
    // RMRS-GTR-371 Table 8a: characteristic SAV 1764 1/ft, relative packing ratio 2.35
    const SurfaceFuelBed fm10 = computeFuelBed(fuelModel10(), FuelMoisture());
    CHECK(std::fabs(fm10.sigma - 1764.0) < 1.0);
    CHECK(std::fabs(fm10.relativePacking - 2.35) < 0.01);
}

TEST(canopy_shelters_the_wind) {
    // No canopy, or a sparse one (crown fill <= 5 %): the unsheltered factor from fuel depth
    CHECK(std::fabs(windAdjustmentFactor(1.0, 0, 0, 0) - unshelteredWaf(1.0)) < 1e-12);
    CHECK(std::fabs(windAdjustmentFactor(1.0, 10, 60, 0.5) - unshelteredWaf(1.0)) < 1e-12);
    // 60 % cover, crown ratio 0.5, 60 ft trees: f = 0.1, WAF = 0.555 / (sqrt(6) ln(41.6 / 7.8)) = 0.135
    const double sheltered = windAdjustmentFactor(1.0, 60, 60, 0.5);
    CHECK(std::fabs(sheltered - 0.135) < 0.002);
    CHECK(sheltered < unshelteredWaf(1.0));
}
