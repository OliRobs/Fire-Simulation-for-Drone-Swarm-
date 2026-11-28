// ============================================================================
// Tests: the run-time grid (grid.h). Every module that is sized by the grid is run here on grids other than
// the app's 1200 x 1000: small ones and non-square ones (in cells and in metres per cell), so the GridSpec
// parameter is really used and nothing still assumes the default size.
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <vector>

#include "boustrophedon.h"
#include "distance_transform.h"
#include "drone.h"
#include "fire_model.h"
#include "grid.h"
#include "mission_controller.h"
#include "mission_planner.h"
#include "perimeter.h"
#include "safe_planner.h"
#include "topography.h"
#include "viewport.h"
#include "voronoi.h"
#include "wind_overlay.h"
#include "wind_solver.h"
#include "test_framework.h"
#include "test_support.h"

namespace {
    // True if calling f throws std::invalid_argument
    template <class F>
    bool throwsInvalidArgument(F f) {
        try { f(); }
        catch (const std::invalid_argument&) { return true; }
        return false;
    }

    // The cells of `big` in the window whose top-left corner is (ox, oy) and whose size is `small`'s, laid out
    // on `small`: so a result on a small grid can be compared with the same part of a result on a big one
    template <class T>
    std::vector<T> crop(const GridSpec& big, const std::vector<T>& values, const GridSpec& small, int ox, int oy) {
        std::vector<T> out(small.count());
        for (int y = 0; y < small.height; ++y)
            for (int x = 0; x < small.width; ++x) out[small.index(x, y)] = values[big.index(ox + x, oy + y)];
        return out;
    }
}

// ----------------------------------------------------------------------------
// GridSpec itself
// ----------------------------------------------------------------------------
TEST(grid_spec_indexes_rows_and_columns) {
    const GridSpec g{ 7, 5, 10.0f, 12.0f };
    CHECK(g.count() == 35);
    CHECK(g.index(0, 0) == 0 && g.index(6, 0) == 6 && g.index(0, 1) == 7 && g.index(6, 4) == 34);
    for (int i = 0; i < 35; ++i) {
        CHECK(g.index(g.xOf(i), g.yOf(i)) == i);
        const Vec2 c = g.centre(i);
        CHECK(c.x == g.xOf(i) + 0.5f && c.y == g.yOf(i) + 0.5f);
    }
    CHECK(g.inside(0, 0) && g.inside(6, 4));
    CHECK(!g.inside(-1, 0) && !g.inside(7, 0) && !g.inside(0, -1) && !g.inside(0, 5));
    const GridRect r = g.rect();
    CHECK(r.x0 == 0 && r.y0 == 0 && r.x1 == 7 && r.y1 == 5);
    CHECK(g == (GridSpec{ 7, 5, 10.0f, 12.0f }) && g != (GridSpec{ 7, 5, 10.0f, 10.0f }));
}

// Public entry points refuse per-cell vectors built for another grid, instead of reading out of bounds
TEST(grid_size_mismatches_are_refused) {
    const GridSpec small{ 120, 100 };
    std::vector<Cell> wrongCells(small.count() + 1);
    std::vector<std::uint8_t> wrongMask(small.count() - 1, 0);
    CHECK(throwsInvalidArgument([&] { FireSimulation sim(small, wrongCells); }));
    CHECK(throwsInvalidArgument([&] { Perimeter p(small, wrongMask); }));
    CHECK(throwsInvalidArgument([&] { DistanceTransform d(small); d.compute(wrongMask); }));
    CHECK(throwsInvalidArgument([&] { computeTopography(small, wrongCells, false); }));
    CHECK(throwsInvalidArgument([&] { solveTerrainWind(small, wrongCells, 20.0f, 270.0f); }));
    // A wind field for another grid
    std::vector<Cell> cells(small.count());
    WindField field;
    field.toward.assign(10, Vec2(1.0f, 0.0f));
    field.speedKmh.assign(10, 20.0f);
    CHECK(throwsInvalidArgument([&] { FireSimulation sim(small, cells, FireWeather(), SIM_MINUTES_PER_TICK, &field); }));
    // The planner's two distance fields on different grids
    DistanceTransform live(small), perim(GridSpec{ 100, 120 });
    SafePlanner planner(live);
    std::mt19937 rng(1);
    CHECK(throwsInvalidArgument([&] { planMission(live, perim, planner, rng); }));
    // The right sizes are accepted
    CHECK(!throwsInvalidArgument([&] { FireSimulation sim(small, cells); }));
}

// ----------------------------------------------------------------------------
// Distance transform
// ----------------------------------------------------------------------------
// The GridSpec constructor gives the same field as the width/height one, and remembers its grid
TEST(distance_transform_on_a_non_square_grid) {
    const GridSpec g{ 150, 60, 10.0f, 25.0f };
    const std::vector<std::uint8_t> fire = discFire(g, 40, 30, 6);
    DistanceTransform a(g, 20), b(g.width, g.height, 20);
    a.compute(fire);
    b.compute(fire);
    CHECK(a.grid() == g);
    bool same = true;
    for (int i = 0; i < static_cast<int>(g.count()); ++i) same = same && a.squaredDistance(i) == b.squaredDistance(i);
    CHECK(same);
    // A cell 10 cells east of the disc's edge is 10 cells away (distances are in cells, whatever the cell size)
    CHECK(std::fabs(std::sqrt(a.squaredDistance(g.index(40 + 6 + 10, 30))) - 10.0) < 1e-9);
}

// ----------------------------------------------------------------------------
// Fire model
// ----------------------------------------------------------------------------
// A fire on a small grid burns exactly the cells the same fire burns on the full grid (as long as it has
// not reached the small grid's edge): the grid size changes nothing but where the map ends
TEST(fire_on_a_small_grid_matches_the_full_grid) {
    const GridSpec small{ 120, 100 };
    const int ox = CX - 60, oy = CY - 50;   // The small grid is the middle of the full one
    std::vector<Cell> bigCells(N), smallCells(small.count());
    FireSimulation big(G, bigCells), little(small, smallCells);
    big.ignite(at(CX, CY));
    little.ignite(small.index(60, 50));
    for (int t = 0; t < 240; ++t) { big.step(); little.step(); }   // 20 minutes: still well inside the window
    const GridRect b = little.fireBounds();
    CHECK(b.x0 > 0 && b.y0 > 0 && b.x1 < small.width && b.y1 < small.height);
    CHECK(crop(G, big.fireMask(), small, ox, oy) == little.fireMask());
    int burned = 0;
    for (std::uint8_t v : little.fireMask()) burned += v;
    std::printf("  (%d cells burned on both grids)\n", burned);
    CHECK(burned > 50);
}

// Non-square cells (10 m east-west, 20 m north-south) on a non-square grid: with no wind and no slope the fire
// covers the same ground distance east and north, so half as many cells north as east
TEST(fire_on_non_square_cells_spreads_evenly_on_the_ground) {
    const GridSpec g{ 200, 120, 10.0f, 20.0f };
    std::vector<Cell> cells(g.count());
    for (Cell& c : cells) c.fuelModel = 109;   // GR9: fast enough with no wind to cover some ground
    FireSimulation sim(g, cells, calm());
    sim.ignite(g.index(100, 60));
    for (int t = 0; t < 600; ++t) sim.step();
    int east = 0, north = 0;
    while (sim.fireMask()[g.index(100 + east + 1, 60)]) ++east;
    while (sim.fireMask()[g.index(100, 60 - north - 1)]) ++north;
    std::printf("  (east %d cells = %d m, north %d cells = %d m)\n", east, east * 10, north, north * 20);
    CHECK(east > 10);
    CHECK(std::abs(east * 10 - north * 20) <= 20);   // Within one 20 m cell
}

// ----------------------------------------------------------------------------
// Topography
// ----------------------------------------------------------------------------
// A plane on a non-square grid of non-square cells: the same slope and uphill direction everywhere
TEST(topography_on_a_non_square_grid) {
    const GridSpec g{ 130, 70, 10.0f, 12.0f };
    std::vector<Cell> cells(g.count());
    // Rises 0.2 m per metre east and 0.1 m per metre south
    for (int y = 0; y < g.height; ++y)
        for (int x = 0; x < g.width; ++x) cells[g.index(x, y)].elevation = 0.2f * x * 10.0f + 0.1f * y * 12.0f;
    computeTopography(g, cells, false);
    const float expected = static_cast<float>(std::atan(std::sqrt(0.2 * 0.2 + 0.1 * 0.1)) * 57.2957795);
    bool ok = true;
    for (int y = 1; y < g.height - 1; ++y)
        for (int x = 1; x < g.width - 1; ++x) {
            const Cell& c = cells[g.index(x, y)];
            ok = ok && std::fabs(c.slopeDeg - expected) < 0.01f && std::fabs(c.upslopeX - 0.894f) < 0.01f
                 && std::fabs(c.upslopeY - 0.447f) < 0.01f;
        }
    CHECK(ok);
}

// ----------------------------------------------------------------------------
// Terrain wind and the wind arrows
// ----------------------------------------------------------------------------
// The wind solver on a small, non-square grid of non-square cells: the field is on that grid, flat ground
// keeps the open wind, and a hill still speeds the wind up over its top
TEST(terrain_wind_on_a_small_non_square_grid) {
    const GridSpec g{ 160, 100, 30.0f, 25.0f };
    WindSolverSettings s = windTestSettings();
    s.columnSizeM = 300.0f;
    std::vector<Cell> flat(g.count());
    const WindField f = solveTerrainWind(g, flat, 20.0f, 225.0f, s);
    CHECK(f.grid == g && f.speedKmh.size() == g.count() && f.toward.size() == g.count());
    CHECK(std::fabs(f.speedKmh[g.index(80, 50)] - 20.0f) < 1e-3f);

    std::vector<Cell> hill(g.count());
    for (int y = 0; y < g.height; ++y)
        for (int x = 0; x < g.width; ++x) {
            const double dx = (x - 80) * 30.0, dy = (y - 50) * 25.0;   // Metres from the summit
            hill[g.index(x, y)].elevation = static_cast<float>(250.0 * std::exp(-(dx * dx + dy * dy) / (2.0 * 700.0 * 700.0)));
        }
    const WindField h = solveTerrainWind(g, hill, 20.0f, 270.0f, s);
    std::printf("  (%d x %d x %d mesh, %d CG iterations, crest %.2f km/h)\n", h.columnsX, h.columnsY, h.layers,
                h.iterations, h.speedKmh[g.index(80, 50)]);
    CHECK(h.converged);
    CHECK(h.speedKmh[g.index(80, 50)] > 1.05f * 20.0f);
}

// Arrows on a small grid stay on it, and the lattice spacing is capped by the small grid's longer side
TEST(wind_arrows_on_a_small_grid) {
    const GridSpec g{ 96, 40 };
    WindField f;
    f.grid = g;
    f.toward.assign(g.count(), windTowardGrid(270.0f));
    f.speedKmh.assign(g.count(), 20.0f);
    const Viewport v(static_cast<float>(g.width * CELL_SIZE), static_cast<float>(g.height * CELL_SIZE));
    const std::vector<WindArrow> arrows = buildWindArrows(f, v, 1.0f, 20.0f, 20.0f);
    CHECK(!arrows.empty());
    bool inside = true;
    for (const WindArrow& a : arrows) {
        const float cx = 0.5f * (a.x0 + a.x1), cy = 0.5f * (a.y0 + a.y1);   // Arrow centre, world units
        inside = inside && cx >= 0 && cx < g.width * CELL_SIZE && cy >= 0 && cy < g.height * CELL_SIZE;
    }
    CHECK(inside);
    // However far out the view is, the lattice never spaces arrows wider than the grid's longer side (128 >= 96)
    CHECK(windArrowSpacing(g, 1000.0f) == 128.0f * CELL_SIZE);
}

// ----------------------------------------------------------------------------
// Planning: perimeter, safe planner, Voronoi, sweep plan, drones
// ----------------------------------------------------------------------------
// The planner's distances are hundreds of cells (the perimeter is 1.2 km out, the search area 2.1 km), so the
// "small" planning grid is 520 x 460 cells: a quarter of the default grid, and not square.
namespace {
    const GridSpec PLAN_GRID{ 520, 460 };
}

// The perimeter on a small grid: same ring and distances as the same fire on the full grid
TEST(perimeter_on_a_small_grid_matches_the_full_grid) {
    const int ox = CX - PLAN_GRID.width / 2, oy = CY - PLAN_GRID.height / 2;
    const std::vector<std::uint8_t> bigFire = discFire(CX, CY, 40);
    const std::vector<std::uint8_t> smallFire = discFire(PLAN_GRID, PLAN_GRID.width / 2, PLAN_GRID.height / 2, 40);
    Perimeter big(G, bigFire), small(PLAN_GRID, smallFire);
    CHECK(small.grid() == PLAN_GRID);
    // Every distance the planner reads (out to SEARCH_DISTANCE) is the same, cell for cell
    constexpr double reach = (SEARCH_DISTANCE + 0.5) * (SEARCH_DISTANCE + 0.5);
    bool same = true;
    for (int y = 0; y < PLAN_GRID.height; ++y)
        for (int x = 0; x < PLAN_GRID.width; ++x) {
            const double a = big.distance().squaredDistance(at(ox + x, oy + y));
            const double b = small.distance().squaredDistance(PLAN_GRID.index(x, y));
            if (a < reach || b < reach) same = same && a == b;
        }
    CHECK(same);
    CHECK(big.ringCells().size() == small.ringCells().size());
    // The traced perimeter line is one loop on both
    CHECK(tracePerimeterLine(small.distance()).size() == 1);
}

// The safe planner on a non-square grid: a leg round a fire avoids it and stays on the grid
TEST(safe_planner_on_a_non_square_grid) {
    const GridSpec g{ 220, 90 };
    const std::vector<std::uint8_t> fire = discFire(g, 110, 45, 20);
    DistanceTransform live(g, LIVE_EDT_RANGE);
    live.compute(fire);
    SafePlanner sp(live);
    const std::vector<Vec2> path = sp.leg(Vec2(10.5f, 45.5f), Vec2(209.5f, 45.5f));
    CHECK(!path.empty());
    Vec2 prev(10.5f, 45.5f);
    bool ok = true;
    for (const Vec2& p : path) {
        ok = ok && !crossesFire(g, fire, prev, p) && p.x >= 0.0f && p.x <= g.width && p.y >= 0.0f && p.y <= g.height;
        prev = p;
    }
    CHECK(ok);
    CHECK(std::fabs(path.back().x - 209.5f) < 1e-4f && std::fabs(path.back().y - 45.5f) < 1e-4f);
    CHECK(!sp.cellOk(g.width, 10) && !sp.cellOk(10, g.height) && sp.cellOk(g.width - 1, 0));
}

// Voronoi on a small non-square grid: every cell gets a region, every region gets cells
TEST(voronoi_on_a_non_square_grid) {
    const GridSpec g{ 90, 40 };
    std::vector<int> cells;
    for (int i = 0; i < static_cast<int>(g.count()); ++i) cells.push_back(i);
    std::mt19937 rng(5);
    const Voronoi v = buildVoronoi(g, cells, 6, rng);
    std::vector<int> count(6, 0);
    bool ok = v.assign.size() == cells.size();
    for (int a : v.assign) { ok = ok && a >= 0 && a < 6; if (a >= 0 && a < 6) ++count[a]; }
    CHECK(ok);
    for (int c : count) CHECK(c > 0);
    for (const Vec2& s : v.sites) CHECK(s.x >= 0.0f && s.x <= g.width && s.y >= 0.0f && s.y <= g.height);
}

// A whole sweep plan on the small planning grid: zones are counted the same as on the full grid, every band
// cell has an owner, and every route stays on the grid
TEST(mission_plan_on_a_small_grid) {
    const GridSpec& g = PLAN_GRID;
    const std::vector<std::uint8_t> fire = discFire(g, g.width / 2, g.height / 2, 30);
    DistanceTransform live(g, LIVE_EDT_RANGE);
    live.compute(fire);
    Perimeter perim(g, fire);
    SafePlanner sp(live);
    std::mt19937 rng(2024);
    const MissionPlan plan = planMission(live, perim.distance(), sp, rng);
    CHECK(plan.zone.size() == g.count() && plan.owner.size() == g.count());

    // The same fire on the full grid classifies the same cells
    const SweepSetup full(discFire(CX, CY, 30));
    const int ox = CX - g.width / 2, oy = CY - g.height / 2;
    CHECK(crop(G, full.plan.zone, g, ox, oy) == plan.zone);

    int band = 0, unowned = 0;
    for (int i = 0; i < static_cast<int>(g.count()); ++i) {
        if (plan.zone[i] != 1) continue;
        ++band;
        if (plan.owner[i] < 0) ++unowned;
    }
    bool onGrid = true;
    for (const auto& route : plan.sweepRoutes)
        for (const Vec2& p : route) onGrid = onGrid && p.x >= 0.0f && p.x <= g.width && p.y >= 0.0f && p.y <= g.height;
    for (const Vec2& p : plan.acoTour) onGrid = onGrid && p.x >= 0.0f && p.x <= g.width && p.y >= 0.0f && p.y <= g.height;
    std::printf("  (%d band cells, %zu lanes, %zu sweep routes, %zu-point ACO tour)\n", band, plan.lanes.size(),
                plan.sweepRoutes.size(), plan.acoTour.size());
    CHECK(band > 0 && unowned == 0);
    CHECK(!plan.lanes.empty() && !plan.sweepRoutes.empty());
    CHECK(onGrid);
    // The ownership map on its own agrees with the planner's
    CHECK(buildOwnerMap(g, plan.lanes, plan.zone, plan.bounds) == plan.owner);
}

// The mission controller flying a growing fire on the small planning grid: no drone crosses the fire, ends a
// frame inside the hard no-fly zone or leaves the grid
TEST(mission_on_a_small_grid_keeps_drones_safe_and_on_the_map) {
    const GridSpec& g = PLAN_GRID;
    std::vector<Cell> cells(g.count());
    FireSimulation sim(g, cells);
    sim.ignite(g.index(g.width / 2, g.height / 2));
    for (int t = 0; t < 300; ++t) sim.step();
    DistanceTransform live(g, LIVE_EDT_RANGE, DistanceTransform::Mode::IncrementalClamped);
    live.compute(sim.fireMask(), sim.fireBounds());
    Perimeter perim(g, sim.fireMask());
    MissionController mission(live, perim.distance(), Vec2(2.5f, g.height - 2.5f));
    int crossings = 0, insideHardZone = 0, offGrid = 0;
    for (int frame = 0; frame < 600; ++frame) {
        sim.step();
        if (sim.takeFireGrew()) {
            live.compute(sim.fireMask(), sim.fireBounds());
            perim.update(sim.fireMask(), sim.fireBounds());
            mission.markFireChanged();
        }
        std::vector<std::pair<int, Vec2>> before;
        for (const Drone& d : mission.drones()) before.emplace_back(d.id, d.pos);
        mission.update();
        for (const Drone& d : mission.drones()) {
            for (const auto& [id, was] : before)
                if (id == d.id && crossesFire(g, sim.fireMask(), was, d.pos)) ++crossings;
            const int x = static_cast<int>(d.pos.x), y = static_cast<int>(d.pos.y);
            if (!g.inside(x, y)) { ++offGrid; continue; }
            if (live.squaredDistance(g.index(x, y)) < HARD_NO_FLY_SQ - 1e-6) ++insideHardZone;
        }
    }
    std::printf("  (%zu drones; moves over fire %d, inside hard zone %d, off the grid %d)\n", mission.drones().size(),
                crossings, insideHardZone, offGrid);
    CHECK(!mission.drones().empty());
    CHECK(crossings == 0 && insideHardZone == 0 && offGrid == 0);
}

// The fire-repulsion safety net keeps a drone on a small grid when it pushes it out towards the grid's edge
TEST(keep_clear_of_fire_stays_on_a_small_grid) {
    const GridSpec g{ 60, 40 };
    // Fire along the whole west edge; a drone just east of it, inside the hard zone
    std::vector<std::uint8_t> fire(g.count(), 0);
    for (int y = 0; y < g.height; ++y) fire[g.index(0, y)] = 1;
    DistanceTransform live(g, LIVE_EDT_RANGE);
    live.compute(fire);
    Drone d;
    d.pos = Vec2(2.5f, 20.5f);
    keepClearOfFire(d, live);
    CHECK(d.pos.x >= 0.5f && d.pos.x <= g.width - 0.5f && d.pos.y >= 0.5f && d.pos.y <= g.height - 0.5f);
    CHECK(live.squaredDistance(g.index(static_cast<int>(d.pos.x), static_cast<int>(d.pos.y))) >= HARD_NO_FLY_SQ);
}

// ----------------------------------------------------------------------------
// Viewport
// ----------------------------------------------------------------------------
// The most zoomed-in view shows VIEW_MIN_CELLS cells across, whatever the map's width
TEST(viewport_max_zoom_follows_the_map_width) {
    for (int cells : { 120, 1200, 2000 }) {
        Viewport v(static_cast<float>(cells * CELL_SIZE), 500.0f);
        v.zoomAt(10.0f, 10.0f, 1e6f);
        CHECK(std::fabs(v.width() - VIEW_MIN_CELLS * CELL_SIZE) < 1e-3f);
    }
}
