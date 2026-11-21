// ============================================================================
// Frame benchmark (target wf_bench): main.cpp's frame loop without SFML, timing each phase.
//
// Run:  ./build/wf_bench [frames] [--expect <checksum>]      (default 2400 frames = 200 simulated minutes)
//
// Synthetic map: GR2 grass everywhere (DEFAULT_FUEL_MODEL), one ignition at the centre of the grid, the
// launch point in the bottom-left corner as in main.cpp. Grid size and frame length come from config.h.
//
// Prints:
//  - per-frame times of each phase (fire step, live distance transform, perimeter update, perimeter ring)
//    and of the mission update (normal frames and re-plans separately);
//  - how often the perimeter moved and the ACO tour changed;
//  - an FNV-1a checksum of the final state. Two builds that give the same checksum gave identical results.
//    With --expect, the program exits with code 1 if the checksum differs (CI uses this; see
//    tools/bench_checksum.txt for the expected value and when to update it).
//
// What the checksum covers, and why distances are clamped: the live distance field is hashed only out to
// LIVE_EDT_RANGE and the perimeter field only out to SEARCH_DISTANCE + 1, the furthest any reader looks.
// Values further out are "very far" for every purpose, so an optimisation that stops storing them exactly
// does not change the checksum, while any change a reader could see does.
// ============================================================================
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "distance_transform.h"
#include "fire_model.h"
#include "mission_controller.h"
#include "perimeter.h"

namespace {
    // ------------------------------------------------------------------------
    // FNV-1a hash over the raw bytes of every value passed to mix(): equal checksums = identical results
    // ------------------------------------------------------------------------
    std::uint64_t g_hash = 1469598103934665603ULL;
    template <class T> void mix(const T& v) {
        const unsigned char* p = reinterpret_cast<const unsigned char*>(&v);
        for (size_t i = 0; i < sizeof(T); ++i) { g_hash ^= p[i]; g_hash *= 1099511628211ULL; }
    }

    // Wall-clock milliseconds between two time points
    using Clock = std::chrono::steady_clock;
    double ms(Clock::time_point a, Clock::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); }

    // Squared distance clamped to (range + 1)^2: everything further than any reader looks hashes the same
    double clamped(double d2, int range) {
        const double cap = (range + 1.0) * (range + 1.0);
        return d2 < cap ? d2 : cap;
    }
}

int main(int argc, char** argv) {
    // --- Arguments: [frames] [--expect <hex checksum>] ---
    int frames = 2400;
    const char* expect = nullptr;
    for (int a = 1; a < argc; ++a) {
        if (std::strcmp(argv[a], "--expect") == 0 && a + 1 < argc) expect = argv[++a];
        else frames = std::atoi(argv[a]);
    }
    if (frames <= 0) { std::fprintf(stderr, "usage: wf_bench [frames] [--expect <checksum>]\n"); return 2; }

    // --- Set-up, as in main.cpp but on a synthetic map (default cells: GR2 grass, flat) ---
    const auto t0 = Clock::now();
    std::vector<Cell> grid(static_cast<size_t>(GRID_WIDTH) * GRID_HEIGHT);
    FireSimulation sim(grid);
    sim.ignite((GRID_HEIGHT / 2) * GRID_WIDTH + GRID_WIDTH / 2);
    DistanceTransform live(GRID_WIDTH, GRID_HEIGHT, LIVE_EDT_RANGE, DistanceTransform::Mode::IncrementalClamped);
    Perimeter perim(sim.fireMask());
    MissionController mission(live, perim.distance(), Vec2(2.5f, GRID_HEIGHT - 2.5f));
    const auto t1 = Clock::now();
    std::printf("grid %d x %d, %.1f s per frame, %d frames (%.0f simulated min), set-up %.0f ms\n", GRID_WIDTH,
                GRID_HEIGHT, SIM_MINUTES_PER_TICK * 60.0f, frames, frames * SIM_MINUTES_PER_TICK, ms(t0, t1));

    // --- Frame loop, in the same order as main.cpp (render calls left out) ---
    double tStep = 0, tLive = 0, tPerim = 0, tRing = 0, tFly = 0, tReplan = 0, worstReplan = 0;
    int grew = 0, ringCalls = 0, replans = 0, perimMoves = 0, tourChanges = 0;
    std::vector<Vec2> lastTour;
    for (int f = 0; f < frames; ++f) {
        auto a = Clock::now();
        sim.step();
        sim.takeDirtyCells();
        auto b = Clock::now();
        tStep += ms(a, b);

        bool ringMoved = false;
        if (sim.takeFireGrew()) {
            ++grew;
            a = Clock::now(); live.compute(sim.fireMask(), sim.fireBounds()); b = Clock::now(); tLive += ms(a, b);
            a = Clock::now(); ringMoved = perim.update(sim.fireMask(), sim.fireBounds()); b = Clock::now(); tPerim += ms(a, b);
            perimMoves += ringMoved;
            mission.markFireChanged();
        }
        // main.cpp re-extracts the ring for drawing when the perimeter moves (and once at the start)
        if (ringMoved || f == 0) {
            a = Clock::now(); const std::vector<int> ring = perim.ringCells(); b = Clock::now();
            tRing += ms(a, b);
            ++ringCalls;
        }

        a = Clock::now(); mission.update(); b = Clock::now();
        const double m = ms(a, b);
        // A re-plan sets the plan-changed flag; every other frame only flies the drones
        if (mission.takePlanChanged()) {
            ++replans;
            tReplan += m;
            if (m > worstReplan) worstReplan = m;
            // Did the ACO tour change with this re-plan? (It should only when the perimeter moved.)
            const std::vector<Vec2>& tour = mission.plan().acoTour;
            bool same = tour.size() == lastTour.size();
            for (size_t k = 0; same && k < tour.size(); ++k) same = tour[k].x == lastTour[k].x && tour[k].y == lastTour[k].y;
            tourChanges += !same;
            lastTour = tour;
        }
        else {
            tFly += m;
        }
    }

    // --- Checksum of the final state ---
    const int N = GRID_WIDTH * GRID_HEIGHT;
    for (int i = 0; i < N; ++i) mix(clamped(live.squaredDistance(i), LIVE_EDT_RANGE));
    for (int i = 0; i < N; ++i) mix(clamped(perim.distance().squaredDistance(i), SEARCH_DISTANCE));
    const MissionPlan& plan = mission.plan();
    for (int i = 0; i < N; ++i) { mix(plan.label[i]); mix(plan.zone[i]); }
    for (int v : plan.owner) mix(v);
    for (int r : perim.ringCells()) mix(r);
    for (const Drone& d : mission.drones()) { mix(d.id); mix(d.pos.x); mix(d.pos.y); }
    for (const auto& route : plan.sweepRoutes) for (const Vec2& p : route) { mix(p.x); mix(p.y); }
    for (const Vec2& p : plan.acoTour) { mix(p.x); mix(p.y); }

    // --- Report ---
    int burned = 0;
    for (std::uint8_t v : sim.fireMask()) burned += v;
    const int flyFrames = frames - replans;
    std::printf("fire grew on %d frames, %d cells burned, %zu drones at the end\n", grew, burned, mission.drones().size());
    std::printf("perimeter moved %d times; ACO tour changed on %d of %d re-plans\n", perimMoves, tourChanges, replans);
    std::printf("per frame [ms]: fire step %.3f, live distance %.3f, perimeter %.3f, ring %.2f (per call, x%d)\n",
                tStep / frames, tLive / frames, tPerim / frames, ringCalls ? tRing / ringCalls : 0.0, ringCalls);
    std::printf("mission [ms]: normal frame %.3f, re-plan average %.1f, worst %.1f\n",
                flyFrames ? tFly / flyFrames : 0.0, replans ? tReplan / replans : 0.0, worstReplan);
    char sum[32];
    std::snprintf(sum, sizeof(sum), "%016llx", static_cast<unsigned long long>(g_hash));
    std::printf("checksum %s\n", sum);

    if (expect && std::strcmp(expect, sum) != 0) {
        std::printf("CHECKSUM MISMATCH: expected %s. Results changed: if that is intended, update "
                    "tools/bench_checksum.txt and say why in the pull request.\n", expect);
        return 1;
    }
    return 0;
}
