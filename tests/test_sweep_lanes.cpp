// ============================================================================
// Tests: Sweep lanes and sectors (boustrophedon.h)
// Part of wf_tests (see test_framework.h for how to run one file or one test).
// ============================================================================
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <unordered_map>
#include <vector>

#include "boustrophedon.h"
#include "distance_transform.h"
#include "mission_planner.h"
#include "perimeter.h"
#include "safe_planner.h"
#include "test_framework.h"
#include "test_support.h"

// ============================================================================
// Sweep lanes and sectors
// ============================================================================
namespace {
    // Distance from p to the nearest run of any lane
    float nearestLane(const std::vector<SweepLane>& lanes, Vec2 p) {
        float best = std::numeric_limits<float>::max();
        for (const SweepLane& l : lanes)
            for (size_t i = 0; i + 1 < l.pts.size(); i += 2) {
                const Vec2 a = l.pts[i], ab = l.pts[i + 1] - a, ap = p - a;
                const float l2 = ab.x * ab.x + ab.y * ab.y;
                const float t = l2 > 1e-9f ? std::clamp((ap.x * ab.x + ap.y * ab.y) / l2, 0.0f, 1.0f) : 0.0f;
                best = std::min(best, std::sqrt(sqDist(p, a + ab * t)));
            }
        return best;
    }

    // Count of inner-band cells that no lane sees (further than LANE_COVER_RADIUS from every lane)
    int unseenInnerCells(const MissionPlan& plan) {
        int unseen = 0;
        for (int y = plan.bounds.y0; y < plan.bounds.y1; ++y)
            for (int x = plan.bounds.x0; x < plan.bounds.x1; ++x)
                if (plan.zone[at(x, y)] == 1 && nearestLane(plan.lanes, Vec2(x + 0.5f, y + 0.5f)) > LANE_COVER_RADIUS + 1e-3f)
                    ++unseen;
        return unseen;
    }

}


// ----------------------------------------------------------------------------
// Original (pre-2026-11-20) gap fill and perimeter-line trace, kept word for word as references.
// The production versions only skip work that cannot change the answer (gap fill: later passes look only at
// the previous pass's patch cells; trace: each sample's distance is worked out once instead of four times), so
// their output must match these exactly, point for point.
// ----------------------------------------------------------------------------
namespace ref {
    constexpr float RAY_STEP = 0.5f;   // Cells between samples along a lane (two per cell, so no cell is skipped)

    // Distance (cells) from the committed fire at cell (x, y); far away off the grid
    float distAt(const DistanceTransform& edt, int x, int y) {
        if (x < 0 || x >= GRID_WIDTH || y < 0 || y >= GRID_HEIGHT) return 1e9f;
        return static_cast<float>(std::sqrt(edt.squaredDistance(y * GRID_WIDTH + x)));
    }

    // ------------------------------------------------------------------------
    // Downhill direction of the distance field at p (towards the nearest fire), unit length.
    // Central differences two cells either side smooth out the cell-to-cell steps of the field. On the ridge
    // between two fire edges (in a bay) the two sides cancel and the gradient is ill-defined; then `fallback`
    // (the inward normal of the perimeter line, or another sensible direction) is used instead.
    // ------------------------------------------------------------------------
    Vec2 downhill(const DistanceTransform& edt, Vec2 p, Vec2 fallback) {
        const int x = static_cast<int>(std::floor(p.x)), y = static_cast<int>(std::floor(p.y));
        const Vec2 g((distAt(edt, x + 2, y) - distAt(edt, x - 2, y)) * 0.25f,
                     (distAt(edt, x, y + 2) - distAt(edt, x, y - 2)) * 0.25f);
        const float gl = std::sqrt(g.x * g.x + g.y * g.y);
        // An exact distance field has gradient length 1 away from ridges; well below that means a ridge
        // (or the edge of the transform's window, where the far side reads "very far")
        if (gl < 0.5f || gl > 2.0f) {
            const float fl = std::sqrt(fallback.x * fallback.x + fallback.y * fallback.y);
            return fl > 1e-6f ? fallback / fl : Vec2(0.0f, 1.0f);
        }
        return g * (-1.0f / gl);
    }

    // ------------------------------------------------------------------------
    // Inner-band runs along the line p + dir * t for t in [t0, t1], as run end points in increasing t.
    // Sampled every RAY_STEP. When stopAtFire is set the walk also stops once it is clearly past the band:
    // closer to the committed fire than the standoff (the band can't continue on the far side of the fire for
    // this lane), or back out beyond the perimeter line after having been inside.
    // ------------------------------------------------------------------------
    std::vector<Vec2> bandRuns(const DistanceTransform& edt, const std::vector<std::uint8_t>& zone,
                               Vec2 p, Vec2 dir, float t0, float t1, bool stopAtFire) {
        std::vector<Vec2> pts;
        bool inRun = false;
        Vec2 last;
        for (float t = t0; t <= t1 + 1e-4f; t += RAY_STEP) {
            const Vec2 q = p + dir * t;
            const int x = static_cast<int>(std::floor(q.x)), y = static_cast<int>(std::floor(q.y));
            if (x < 0 || x >= GRID_WIDTH || y < 0 || y >= GRID_HEIGHT) {
                if (t > 0.0f) break;   // Ran off the map
                continue;
            }
            if (stopAtFire && t > 0.0f) {
                const float d = distAt(edt, x, y);
                if (d < MIN_FIRE_DISTANCE - 1.0f) break;                    // Past the standoff, at the fire
                if (t > 2.0f && d > PERIMETER_DISTANCE + 1.5f) break;       // Left the band outwards again
            }
            const bool in = zone[y * GRID_WIDTH + x] == 1;
            if (in && !inRun) { pts.push_back(q); inRun = true; }           // A run starts
            if (!in && inRun) { pts.push_back(last); inRun = false; }       // A run ends at the previous sample
            if (in) last = q;
        }
        if (inRun) pts.push_back(last);
        return pts;
    }

    // Distance from point p to segment a-b
    float segDist(Vec2 p, Vec2 a, Vec2 b) {
        const Vec2 ab = b - a, ap = p - a;
        const float l2 = ab.x * ab.x + ab.y * ab.y;
        const float t = l2 > 1e-9f ? std::clamp((ap.x * ab.x + ap.y * ab.y) / l2, 0.0f, 1.0f) : 0.0f;
        return std::sqrt(sqDist(p, a + ab * t));
    }

    // Distance from p to the nearest point of a lane (its runs only, not the gaps between them)
    float laneDist(const SweepLane& lane, Vec2 p) {
        float best = std::numeric_limits<float>::max();
        for (size_t i = 0; i + 1 < lane.pts.size(); i += 2) best = std::min(best, segDist(p, lane.pts[i], lane.pts[i + 1]));
        return best;
    }

    // Marks every cell (inside area, as a local bitmap) whose centre is within LANE_COVER_RADIUS of the lane's
    // runs. Each run is sampled every RAY_STEP and a disc stamped at each sample; a cell within the radius of a
    // sample is within the radius of the run, so a stamped cell is always really covered.
    void stampLane(const SweepLane& lane, const GridRect& area, std::vector<std::uint8_t>& covered) {
        const float r = LANE_COVER_RADIUS, r2 = r * r;
        const int w = area.x1 - area.x0;
        for (size_t i = 0; i + 1 < lane.pts.size(); i += 2) {
            const Vec2 a = lane.pts[i], b = lane.pts[i + 1];
            const float len = std::sqrt(sqDist(a, b));
            const int steps = std::max(1, static_cast<int>(std::ceil(len / RAY_STEP)));
            for (int s = 0; s <= steps; ++s) {
                const Vec2 c = a + (b - a) * (static_cast<float>(s) / steps);
                const int x0 = std::max(area.x0, static_cast<int>(std::floor(c.x - r)));
                const int x1 = std::min(area.x1 - 1, static_cast<int>(std::floor(c.x + r)));
                const int y0 = std::max(area.y0, static_cast<int>(std::floor(c.y - r)));
                const int y1 = std::min(area.y1 - 1, static_cast<int>(std::floor(c.y + r)));
                for (int y = y0; y <= y1; ++y)
                    for (int x = x0; x <= x1; ++x)
                        if (sqDist(Vec2(x + 0.5f, y + 0.5f), c) <= r2) covered[(y - area.y0) * w + (x - area.x0)] = 1;
            }
        }
    }

int fillCoverageGapsOriginal(std::vector<SweepLane>& lanes, const DistanceTransform& perimEdt,
                     const std::vector<std::uint8_t>& zone, const GridRect& area) {
    if (area.empty() || lanes.empty()) return 0;
    const int w = area.x1 - area.x0, h = area.y1 - area.y0;
    std::vector<std::uint8_t> covered(static_cast<size_t>(w) * h, 0);
    for (const SweepLane& l : lanes) stampLane(l, area, covered);

    constexpr int MAX_PASSES = 6;
    int added = 0;
    for (int pass = 0; pass < MAX_PASSES; ++pass) {
        // --- Find patches of unseen inner cells (8-connected flood fill) ---
        std::vector<int> patchOf(static_cast<size_t>(w) * h, -1);
        std::vector<std::vector<int>> patches;   // Local indices of the cells in each patch
        for (int i = 0; i < w * h; ++i) {
            const int gx = area.x0 + i % w, gy = area.y0 + i / w;
            if (covered[i] || patchOf[i] >= 0 || zone[gy * GRID_WIDTH + gx] != 1) continue;
            const int id = static_cast<int>(patches.size());
            patches.emplace_back();
            std::vector<int> stack{ i };
            patchOf[i] = id;
            while (!stack.empty()) {
                const int c = stack.back();
                stack.pop_back();
                patches[id].push_back(c);
                const int cx = c % w, cy = c / w;
                for (int k = 0; k < NUM_DIRS; ++k) {
                    const int nx = cx + DX[k], ny = cy + DY[k];
                    if (nx < 0 || nx >= w || ny < 0 || ny >= h) continue;
                    const int ni = ny * w + nx;
                    if (covered[ni] || patchOf[ni] >= 0) continue;
                    if (zone[(area.y0 + ny) * GRID_WIDTH + area.x0 + nx] != 1) continue;
                    patchOf[ni] = id;
                    stack.push_back(ni);
                }
            }
        }
        if (patches.empty()) break;   // Everything is seen

        // --- One short lane through each patch ---
        const size_t mainCount = lanes.size();
        std::vector<SweepLane> extra;
        for (const std::vector<int>& patch : patches) {
            // Patch centre: the patch cell nearest the mean of its cells (the mean itself can lie outside a
            // curved patch)
            Vec2 mean(0.0f, 0.0f);
            for (int c : patch) mean += Vec2(area.x0 + c % w + 0.5f, area.y0 + c / w + 0.5f);
            mean = mean / static_cast<float>(patch.size());
            Vec2 centre = mean;
            float bd = std::numeric_limits<float>::max();
            for (int c : patch) {
                const Vec2 p(area.x0 + c % w + 0.5f, area.y0 + c / w + 0.5f);
                if (sqDist(p, mean) < bd) { bd = sqDist(p, mean); centre = p; }
            }

            // Same direction as the main lanes: down the distance-field gradient at the centre. Fallback on a
            // ridge: towards the nearest main lane's inner end direction (any lane direction will do there).
            const SweepLane& nearLane = *std::min_element(lanes.begin(), lanes.begin() + mainCount,
                [&](const SweepLane& a, const SweepLane& b) { return laneDist(a, centre) < laneDist(b, centre); });
            const Vec2 laneDir = nearLane.pts.back() - nearLane.pts.front();
            const Vec2 dir = downhill(perimEdt, centre, laneDir.x * laneDir.x + laneDir.y * laneDir.y > 1e-6f
                                                         ? laneDir : Vec2(0.0f, 1.0f));

            // How far along that line the patch reaches (cells close enough to the line to be seen from it)
            float tMin = 0.0f, tMax = 0.0f;
            for (int c : patch) {
                const Vec2 rel = Vec2(area.x0 + c % w + 0.5f, area.y0 + c / w + 0.5f) - centre;
                const float t = rel.x * dir.x + rel.y * dir.y;
                const float perp = std::fabs(rel.x * dir.y - rel.y * dir.x);
                if (perp > LANE_COVER_RADIUS) continue;
                tMin = std::min(tMin, t);
                tMax = std::max(tMax, t);
            }

            SweepLane lane;
            lane.pts = bandRuns(perimEdt, zone, centre, dir, tMin - RAY_STEP, tMax + RAY_STEP, false);
            if (lane.pts.empty()) lane.pts = { centre, centre };   // Centre is an inner cell, so can't happen
            lane.seed = lane.pts.front();
            lane.loop = nearLane.loop;
            lane.gapFill = true;
            extra.push_back(std::move(lane));
        }

        // --- Insert each new lane after the main lane it sits closest to (and after earlier gap lanes there) ---
        for (SweepLane& g : extra) {
            stampLane(g, area, covered);
            const Vec2 mid = (g.pts.front() + g.pts.back()) * 0.5f;
            size_t parent = 0;
            float pd = std::numeric_limits<float>::max();
            for (size_t i = 0; i < lanes.size(); ++i) {
                if (lanes[i].gapFill || lanes[i].loop != g.loop) continue;
                const float d = laneDist(lanes[i], mid);
                if (d < pd) { pd = d; parent = i; }
            }
            size_t at = parent + 1;
            while (at < lanes.size() && lanes[at].gapFill) ++at;
            lanes.insert(lanes.begin() + static_cast<std::ptrdiff_t>(at), std::move(g));
            ++added;
        }
    }
    return added;
}

    // Corner offsets of a square, in walking order: top-left, top-right, bottom-right, bottom-left
    constexpr int SQ_DX[4] = { 0, 1, 1, 0 };
    constexpr int SQ_DY[4] = { 0, 0, 1, 1 };

    // Id of the square side from corner c to corner c+1 (mod 4) of the square whose top-left sample is (x, y).
    // Horizontal sides and vertical sides get even and odd ids, keyed by the side's top / left sample, so the
    // side shared by two squares gets the same id from both. Samples run from -1 to GRID_WIDTH (one ring of
    // padding round the grid), hence the +1 offsets.
    std::int64_t sideId(int x, int y, int c) {
        const std::int64_t W = GRID_WIDTH + 2;
        auto key = [&](int sx, int sy, int vertical) {
            return ((static_cast<std::int64_t>(sy) + 1) * W + (sx + 1)) * 2 + vertical;
        };
        switch (c) {
        case 0:  return key(x, y, 0);         // Top side
        case 1:  return key(x + 1, y, 1);     // Right side
        case 2:  return key(x, y + 1, 0);     // Bottom side
        default: return key(x, y, 1);         // Left side
        }
    }

    // One traced piece of line inside a square, from one side crossing to another
    struct Piece {
        Vec2 from, to;
        std::int64_t toSide;   // Side the piece ends on, i.e. where the next piece starts
        bool used = false;
    };

std::vector<std::vector<Vec2>> tracePerimeterLineOriginal(const DistanceTransform& perimEdt, float level) {
    // Only the transform's window holds finite distances; the level line lies well inside it, so squares
    // from one sample outside the window to its far edge are enough
    const GridRect win = perimEdt.window();
    std::vector<std::vector<Vec2>> loops;
    if (win.empty()) return loops;

    // Level-line function at sample (x, y): negative on the fire side of the line (nearer than `level`).
    // Samples off the grid count as far from the fire, so a loop that would leave the map closes along its edge.
    auto f = [&](int x, int y) -> double {
        if (x < 0 || x >= GRID_WIDTH || y < 0 || y >= GRID_HEIGHT) return 1e9;
        return std::sqrt(perimEdt.squaredDistance(y * GRID_WIDTH + x)) - level;
    };

    // Pieces keyed by the side they start on (each side crossing starts exactly one piece)
    std::unordered_map<std::int64_t, Piece> pieces;
    std::vector<std::int64_t> order;   // Start sides in scan order, so tracing is deterministic

    for (int y = win.y0 - 1; y < win.y1; ++y) {
        for (int x = win.x0 - 1; x < win.x1; ++x) {
            double v[4];
            bool in[4];
            int nIn = 0;
            for (int c = 0; c < 4; ++c) {
                v[c] = f(x + SQ_DX[c], y + SQ_DY[c]);
                in[c] = v[c] < 0.0;
                nIn += in[c];
            }
            if (nIn == 0 || nIn == 4) continue;   // The line doesn't pass through this square

            // Crossing point on side c (corner c -> c+1), and whether the walk leaves the fire side there
            Vec2 at[4];
            bool leaving[4], crosses[4];
            for (int c = 0; c < 4; ++c) {
                const int d = (c + 1) % 4;
                crosses[c] = in[c] != in[d];
                leaving[c] = in[c] && !in[d];
                if (!crosses[c]) continue;
                const double t = v[c] / (v[c] - v[d]);   // Where the linear interpolation reaches 0
                const Vec2 a(x + SQ_DX[c] + 0.5f, y + SQ_DY[c] + 0.5f), b(x + SQ_DX[d] + 0.5f, y + SQ_DY[d] + 0.5f);
                at[c] = a + (b - a) * static_cast<float>(t);
            }

            // Pair each leaving crossing with an entering one. With two crossings there is only one way. With
            // four (a saddle: opposite corners on the same side) the centre value decides: if the centre is on
            // the fire side, the two fire corners are joined, so each leaving crossing pairs with the NEXT
            // entering crossing round the square (cutting off the far corners); otherwise with the previous one.
            const bool saddle = nIn == 2 && in[0] == in[2];
            const bool centreIn = (v[0] + v[1] + v[2] + v[3]) * 0.25 < 0.0;
            for (int c = 0; c < 4; ++c) {
                if (!crosses[c] || !leaving[c]) continue;
                int e = -1;
                if (saddle && !centreIn) {
                    for (int k = 1; k < 4 && e < 0; ++k) { const int j = (c + 4 - k) % 4; if (crosses[j] && !leaving[j]) e = j; }
                }
                else {
                    for (int k = 1; k < 4 && e < 0; ++k) { const int j = (c + k) % 4; if (crosses[j] && !leaving[j]) e = j; }
                }
                if (e < 0) continue;   // Cannot happen with a consistent sign pattern
                const std::int64_t start = sideId(x, y, c);
                pieces[start] = Piece{ at[c], at[e], sideId(x, y, e) };
                order.push_back(start);
            }
        }
    }

    // Chain pieces into loops: each piece ends on the side where the next one starts
    for (std::int64_t first : order) {
        auto it = pieces.find(first);
        if (it == pieces.end() || it->second.used) continue;
        std::vector<Vec2> loop;
        while (it != pieces.end() && !it->second.used) {
            it->second.used = true;
            loop.push_back(it->second.from);
            it = pieces.find(it->second.toSide);
        }
        if (loop.size() >= 3) loops.push_back(std::move(loop));
    }

    // Make every loop run anticlockwise on screen (negative area in y-down grid coordinates). The piece rule
    // gives every loop the same sense round the fire, so if the first enclosing loop is the wrong way round,
    // they all are: test the loop with the largest area (always one that encloses fire) and flip all if needed.
    if (!loops.empty()) {
        size_t biggest = 0;
        for (size_t i = 1; i < loops.size(); ++i) {
            if (std::fabs(loopSignedArea(loops[i])) > std::fabs(loopSignedArea(loops[biggest]))) biggest = i;
        }
        if (loopSignedArea(loops[biggest]) > 0.0) {
            for (auto& l : loops) std::reverse(l.begin(), l.end());
        }
    }

    // Start each loop at its topmost point (then leftmost), a stable place for lane numbering to begin
    for (auto& l : loops) {
        const auto top = std::min_element(l.begin(), l.end(), [](Vec2 a, Vec2 b) {
            return a.y < b.y || (a.y == b.y && a.x < b.x);
        });
        std::rotate(l.begin(), top, l.end());
    }
    return loops;
}
}   // namespace ref


// Lane coverage is now worked out per cell (only the samples that can reach each cell are tested). It must mark
// exactly the cells the original disc-per-sample stamp (ref::stampLane, kept above) marks. Random lanes at every
// angle and length, with several runs, axis-aligned and zero-length runs, clipped by the area's edges.
TEST(lane_stamp_matches_disc_stamps) {
    std::mt19937 rng(11);
    std::uniform_real_distribution<float> pos(300.0f, 700.0f), len(0.0f, 300.0f), ang(0.0f, 6.2831853f), gap(0.0f, 40.0f);
    const GridRect area{ 380, 330, 760, 640 };   // Smaller than the spread of the lanes, so many are clipped
    const size_t cells = static_cast<size_t>(area.x1 - area.x0) * (area.y1 - area.y0);
    int lanes = 0, mismatched = 0;
    long long marked = 0;
    for (int k = 0; k < 2000; ++k) {
        SweepLane lane;
        Vec2 p(pos(rng), pos(rng));
        // Direction: every 10th lane axis-aligned (a special case of the row interval), the rest at random
        const float th = k % 10 == 0 ? 1.5707963f * (k / 10 % 4) : ang(rng);
        const Vec2 dir(std::cos(th), std::sin(th));
        const int runs = 1 + k % 3;
        for (int r = 0; r < runs; ++r) {
            // Every 7th run has zero length, every 13th a tiny one (just above the zero-length cut-off)
            const float l = k % 7 == 0 ? 0.0f : (k % 13 == 0 ? 2e-3f : len(rng) / runs);
            lane.pts.push_back(p);
            p = p + dir * l;
            lane.pts.push_back(p);
            p = p + dir * gap(rng);
        }
        std::vector<std::uint8_t> fast(cells, 0), original(cells, 0);
        stampLaneCoverage(lane, area, fast);
        ref::stampLane(lane, area, original);
        mismatched += fast != original;
        for (std::uint8_t v : original) marked += v;
        ++lanes;
    }
    std::printf("  (%d random lanes, %lld cells marked, %d lanes differ)\n", lanes, marked, mismatched);
    CHECK(mismatched == 0);
    CHECK(marked > 0);
}

// The stamp now tests a cell only against the samples within h = sqrt((r + 0.1)^2 - across^2) of its projection
// onto the run (2026-11-23). The cells where skipping a sample could go wrong are the ones right on the swath's
// edge, where the nearest sample is almost exactly LANE_COVER_RADIUS away. Lanes on whole and half-cell
// coordinates put cell centres EXACTLY r = 6.5 cells across from the run (and exactly on samples along it), and
// lanes nudged by a few float ulps put them a hair inside or outside; long lanes far out on a 2000-cell grid
// have the largest rounding. Every stamp must match the disc stamps bit for bit.
SLOW_TEST(lane_stamp_matches_disc_stamps_on_the_edge) {
    std::mt19937 rng(23);
    std::uniform_int_distribution<int> cell(20, 1950);
    std::uniform_int_distribution<int> lenCells(0, 400), nudge(-3, 3), half(0, 1), axis(0, 5);
    const GridRect area{ 0, 0, 2000, 2000 };   // Whole 2000 x 2000 grid, coordinates up to 2000 (most rounding)
    const size_t cells = static_cast<size_t>(area.x1 - area.x0) * (area.y1 - area.y0);
    std::vector<std::uint8_t> fast(cells), original(cells);
    int lanes = 0, mismatched = 0;
    long long marked = 0;
    // Moves a float by k units in the last place (k may be negative), to sit just either side of the edge
    auto ulps = [](float v, int k) {
        for (; k > 0; --k) v = std::nextafter(v, 1e9f);
        for (; k < 0; ++k) v = std::nextafter(v, -1e9f);
        return v;
    };
    for (int k = 0; k < 600; ++k) {
        // Start on a whole or half cell, nudged a few ulps; run along an axis, a diagonal or a 1:2 slope
        const float sx = ulps(cell(rng) + 0.5f * half(rng), nudge(rng));
        const float sy = ulps(cell(rng) + 0.5f * half(rng), nudge(rng));
        static const int DIRS[6][2] = { {1, 0}, {0, 1}, {-1, 0}, {1, 1}, {1, -2}, {2, 1} };
        const int* d = DIRS[axis(rng)];
        const float l = static_cast<float>(lenCells(rng));
        const float norm = std::sqrt(static_cast<float>(d[0] * d[0] + d[1] * d[1]));
        const Vec2 a(sx, sy);
        Vec2 b(std::clamp(sx + d[0] * l / norm, 1.0f, 1999.0f), std::clamp(sy + d[1] * l / norm, 1.0f, 1999.0f));
        b = Vec2(ulps(b.x, nudge(rng)), ulps(b.y, nudge(rng)));
        SweepLane lane;
        lane.pts = { a, b };
        std::fill(fast.begin(), fast.end(), 0);
        std::fill(original.begin(), original.end(), 0);
        stampLaneCoverage(lane, area, fast);
        ref::stampLane(lane, area, original);
        mismatched += fast != original;
        for (std::uint8_t v : original) marked += v;
        ++lanes;
    }
    std::printf("  (%d edge lanes, %lld cells marked, %d lanes differ)\n", lanes, marked, mismatched);
    CHECK(mismatched == 0);
    CHECK(marked > 0);
}

// Gap fill now only re-floods the previous pass's patches (lanes only remove gaps, never add them). The lanes
// it adds must be exactly the ones the original whole-window scan added, in the same places and order.
TEST(gap_fill_matches_full_window_scan) {
    int totalAdded = 0;
    for (const auto& fire : exactnessFires()) {
        SweepSetup s(fire);
        // Rebuild the gap fill's input the way planMission does, then run both versions on copies of it
        const std::vector<SweepLane> lanes = layLanes(tracePerimeterLine(s.perim.distance()), s.perim.distance(),
                                                      s.plan.zone);
        std::vector<SweepLane> fast = lanes, original = lanes;
        const int addedFast = fillCoverageGaps(fast, s.perim.distance(), s.plan.zone, s.plan.bounds);
        const int addedOriginal = ref::fillCoverageGapsOriginal(original, s.perim.distance(), s.plan.zone, s.plan.bounds);
        CHECK(addedFast == addedOriginal);
        CHECK(sameLanes(fast, original));
        totalAdded += addedOriginal;
    }
    std::printf("  (%d gap-fill lanes added over the test fires, identical in both versions)\n", totalAdded);
    CHECK(totalAdded > 0);   // The fires really exercise the gap fill
}

// Original (pre-2026-11-22) ownership flood, kept word for word: seed every cell a lane run passes over, then
// breadth-first through the band, bounds-testing every neighbour. The last stage of buildOwnerMap (band cells cut
// off from every lane go to the nearest lane) is unchanged and left out here, so those cells stay -1.
namespace ref {
    void ownerSeedAndFlood(const std::vector<SweepLane>& lanes, const std::vector<std::uint8_t>& zone,
                           const GridRect& area, std::vector<int>& owner) {
        constexpr float RAY_STEP = 0.5f;
        std::vector<int> queue;
        for (const SweepLane& l : lanes) {
            for (size_t i = 0; i + 1 < l.pts.size(); i += 2) {
                const Vec2 a = l.pts[i], b = l.pts[i + 1];
                const int steps = std::max(1, static_cast<int>(std::ceil(std::sqrt(sqDist(a, b)) / RAY_STEP)));
                for (int s = 0; s <= steps; ++s) {
                    const Vec2 p = a + (b - a) * (static_cast<float>(s) / steps);
                    const int x = static_cast<int>(std::floor(p.x)), y = static_cast<int>(std::floor(p.y));
                    if (x < area.x0 || x >= area.x1 || y < area.y0 || y >= area.y1) continue;
                    const int idx = y * GRID_WIDTH + x;
                    if (zone[idx] != 1 || owner[idx] >= 0) continue;
                    owner[idx] = l.sector;
                    queue.push_back(idx);
                }
            }
        }
        for (size_t head = 0; head < queue.size(); ++head) {
            const int idx = queue[head];
            const int x = idx % GRID_WIDTH, y = idx / GRID_WIDTH;
            for (int k = 0; k < NUM_DIRS; ++k) {
                const int nx = x + DX[k], ny = y + DY[k];
                if (nx < area.x0 || nx >= area.x1 || ny < area.y0 || ny >= area.y1) continue;
                const int ni = ny * GRID_WIDTH + nx;
                if (zone[ni] != 1 || owner[ni] >= 0) continue;
                owner[ni] = owner[idx];
                queue.push_back(ni);
            }
        }
    }
}   // namespace ref

// buildOwnerMap's breadth-first flood now skips the bounds tests for cells inside the area and reserves its queue.
// Every cell the flood reaches must get exactly the sector the original flood gave it (same lane reaches it first),
// and the cells it cannot reach must still be given a sector by the last stage. Run on the test fires with the
// planner's own window and with smaller areas cut out of it, so the area's edge runs through the band (many edge
// cells), plus an island of band cells no lane reaches.
TEST(owner_map_matches_original_flood) {
    long long compared = 0, cutOff = 0, differ = 0, unowned = 0;
    for (const auto& fire : exactnessFires()) {
        SweepSetup s(fire);
        const GridRect w = s.plan.bounds;
        const int mx = (w.x0 + w.x1) / 2, my = (w.y0 + w.y1) / 2;
        // An island of band cells in the window's corner, well away from the band and every lane, so the flood
        // cannot reach it and the last stage has to place it
        std::vector<std::uint8_t> zone = s.plan.zone;
        for (int y = w.y0 + 2; y < w.y0 + 12; ++y)
            for (int x = w.x0 + 2; x < w.x0 + 12; ++x) zone[at(x, y)] = 1;
        for (const GridRect area : { w, GridRect{ w.x0, w.y0, mx, my }, GridRect{ mx - 37, w.y0, w.x1, my + 23 },
                                     GridRect{ w.x0 + 5, my - 60, w.x1 - 5, my + 60 } }) {
            std::vector<int> fast(N, -1), original(N, -1);
            buildOwnerMap(s.plan.lanes, zone, area, fast);
            ref::ownerSeedAndFlood(s.plan.lanes, zone, area, original);
            for (int y = area.y0; y < area.y1; ++y) {
                for (int x = area.x0; x < area.x1; ++x) {
                    const int i = at(x, y);
                    if (zone[i] != 1) { differ += fast[i] != -1; continue; }
                    if (original[i] >= 0) { differ += fast[i] != original[i]; ++compared; }
                    else { ++cutOff; unowned += fast[i] < 0 && !s.plan.lanes.empty(); }
                }
            }
            // Nothing outside the area is written
            for (int i = 0; i < N; ++i) {
                const int x = i % GRID_WIDTH, y = i / GRID_WIDTH;
                if (x < area.x0 || x >= area.x1 || y < area.y0 || y >= area.y1) differ += fast[i] != -1;
            }
        }
    }
    std::printf("  (%lld flooded cells compared, %lld cut off from every lane; %lld differ, %lld left unowned)\n",
                compared, cutOff, differ, unowned);
    CHECK(differ == 0);
    CHECK(unowned == 0);
    CHECK(compared > 100000);
    CHECK(cutOff > 0);   // The cut-down areas really leave band pieces the flood cannot reach
}

// buildOwnerMap can be told how many band cells the area holds (the planner counts them while classifying) and
// then skips its final scan for cut-off cells when the flood has already reached them all. The map must be
// identical with and without the count, both when nothing is cut off (scan skipped) and when an island of band
// cells is cut off from every lane (scan still runs).
TEST(owner_map_band_count_skip_is_exact) {
    int skippedCases = 0, scannedCases = 0;
    long long differ = 0;
    for (const auto& fire : exactnessFires()) {
        SweepSetup s(fire);
        const GridRect w = s.plan.bounds;
        std::vector<std::uint8_t> island = s.plan.zone;   // As owner_map_matches_original_flood: a far island
        for (int y = w.y0 + 2; y < w.y0 + 12; ++y)
            for (int x = w.x0 + 2; x < w.x0 + 12; ++x) island[at(x, y)] = 1;
        for (const std::vector<std::uint8_t>* zone : { &s.plan.zone, &island }) {
            long long band = 0;
            for (int y = w.y0; y < w.y1; ++y)
                for (int x = w.x0; x < w.x1; ++x) band += (*zone)[at(x, y)] == 1;
            const std::vector<int> plain = buildOwnerMap(s.plan.lanes, *zone, w);          // Always scans
            const std::vector<int> counted = buildOwnerMap(s.plan.lanes, *zone, w, band);  // May skip the scan
            for (int i = 0; i < N; ++i) differ += plain[i] != counted[i];
            // Every band cell has an owner either way (so a skipped scan really had nothing left to do)
            long long unowned = 0;
            for (int y = w.y0; y < w.y1; ++y)
                for (int x = w.x0; x < w.x1; ++x) unowned += (*zone)[at(x, y)] == 1 && counted[at(x, y)] < 0;
            CHECK(unowned == 0 || s.plan.lanes.empty());
            (zone == &island ? scannedCases : skippedCases) += 1;
        }
    }
    std::printf("  (%d plain and %d cut-off cases; %lld cells differ)\n", skippedCases, scannedCases, differ);
    CHECK(differ == 0);
    CHECK(skippedCases > 0 && scannedCases > 0);
}

// The planner now builds the ownership map with buildOwnerMapUsingZone, which marks claimed band cells in the
// zone map itself (one byte read per neighbour instead of the zone and owner arrays) and puts the marks back.
// Checked on the test fires, with the planner's window and cut-down areas (many edge cells, pieces cut off from
// every lane) and a far island of band cells:
//  - every cell the flood reaches gets exactly the sector the original flood (ref::ownerSeedAndFlood) gave it;
//  - cut-off band cells still get an owner from the final scan, and nothing outside the area is written;
//  - the zone map is exactly as it was on entry, with and without the band count (which skips the final scan).
TEST(owner_map_using_zone_matches_and_restores_zone) {
    long long compared = 0, differ = 0, unowned = 0, zoneChanged = 0;
    int cases = 0;
    for (const auto& fire : exactnessFires()) {
        SweepSetup s(fire);
        const GridRect w = s.plan.bounds;
        const int mx = (w.x0 + w.x1) / 2, my = (w.y0 + w.y1) / 2;
        std::vector<std::uint8_t> island = s.plan.zone;   // As owner_map_matches_original_flood: a far island
        for (int y = w.y0 + 2; y < w.y0 + 12; ++y)
            for (int x = w.x0 + 2; x < w.x0 + 12; ++x) island[at(x, y)] = 1;
        for (const std::vector<std::uint8_t>* zoneIn : { &s.plan.zone, &island }) {
            for (const GridRect area : { w, GridRect{ w.x0, w.y0, mx, my }, GridRect{ mx - 37, w.y0, w.x1, my + 23 } }) {
                // Band cells in the area, for the form that may skip the final scan
                long long band = 0;
                for (int y = area.y0; y < area.y1; ++y)
                    for (int x = area.x0; x < area.x1; ++x) band += (*zoneIn)[at(x, y)] == 1;
                std::vector<int> original(N, -1);
                ref::ownerSeedAndFlood(s.plan.lanes, *zoneIn, area, original);
                for (const long long count : { -1LL, band }) {
                    std::vector<std::uint8_t> zone = *zoneIn;   // Marked and restored by the call
                    std::vector<int> fast(N, -1);
                    buildOwnerMapUsingZone(s.plan.lanes, zone, area, fast, count);
                    zoneChanged += zone != *zoneIn;
                    ++cases;
                    for (int i = 0; i < N; ++i) {
                        const int x = i % GRID_WIDTH, y = i / GRID_WIDTH;
                        const bool inside = x >= area.x0 && x < area.x1 && y >= area.y0 && y < area.y1;
                        if (!inside || (*zoneIn)[i] != 1) { differ += fast[i] != -1; continue; }
                        if (original[i] >= 0) { differ += fast[i] != original[i]; ++compared; }
                        else unowned += fast[i] < 0 && !s.plan.lanes.empty();
                    }
                }
            }
        }
    }
    std::printf("  (%d cases, %lld flooded cells compared; %lld differ, %lld left unowned, zone changed in %lld)\n",
                cases, compared, differ, unowned, zoneChanged);
    CHECK(differ == 0);
    CHECK(unowned == 0);
    CHECK(zoneChanged == 0);
    CHECK(compared > 100000);
}

// Original (pre-2026-11-25) sector matching, kept word for word: it counted shared cells by scanning every cell of
// the plan window and searched the fleet's ids for each one. The planner now counts over its list of band cells
// with a sector-to-drone table worked out once.
namespace ref {
    std::vector<int> matchSectorsOriginal(const MissionPlan& plan, const std::vector<std::vector<Vec2>>& routes,
                                          const SweepFleet& fleet, const GridRect& area) {
        const int S = static_cast<int>(routes.size()), D = static_cast<int>(fleet.ids.size());
        std::vector<int> result(S, -1);
        std::vector<char> droneUsed(D, 0);
        struct Pair { double score; int sector, drone; };
        std::vector<Pair> pairs;
        const MissionPlan* prev = fleet.previous;
        if (prev && !prev->owner.empty() && !prev->sectorDrone.empty()) {
            std::vector<long long> shared(static_cast<size_t>(S) * D, 0);
            for (int y = area.y0; y < area.y1; ++y) {
                for (int i = y * GRID_WIDTH + area.x0, end = y * GRID_WIDTH + area.x1; i < end; ++i) {
                    const int s = plan.owner[i], ps = prev->owner[i];
                    if (s < 0 || ps < 0 || ps >= static_cast<int>(prev->sectorDrone.size())) continue;
                    const int id = prev->sectorDrone[ps];
                    for (int d = 0; d < D; ++d) if (fleet.ids[d] == id) { ++shared[static_cast<size_t>(s) * D + d]; break; }
                }
            }
            for (int s = 0; s < S; ++s)
                for (int d = 0; d < D; ++d)
                    if (shared[static_cast<size_t>(s) * D + d] > 0) pairs.push_back({ static_cast<double>(shared[static_cast<size_t>(s) * D + d]), s, d });
            std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) {
                return a.score > b.score || (a.score == b.score && (a.sector < b.sector || (a.sector == b.sector && a.drone < b.drone)));
            });
            for (const Pair& p : pairs) {
                if (result[p.sector] >= 0 || droneUsed[p.drone]) continue;
                result[p.sector] = fleet.ids[p.drone];
                droneUsed[p.drone] = 1;
            }
        }
        pairs.clear();
        for (int s = 0; s < S; ++s) {
            if (result[s] >= 0 || routes[s].empty()) continue;
            for (int d = 0; d < D; ++d) {
                if (droneUsed[d]) continue;
                float best = std::numeric_limits<float>::max();
                for (const Vec2& p : routes[s]) best = std::min(best, sqDist(p, fleet.positions[d]));
                pairs.push_back({ best, s, d });
            }
        }
        std::sort(pairs.begin(), pairs.end(), [](const Pair& a, const Pair& b) {
            return a.score < b.score || (a.score == b.score && (a.sector < b.sector || (a.sector == b.sector && a.drone < b.drone)));
        });
        for (const Pair& p : pairs) {
            if (result[p.sector] >= 0 || droneUsed[p.drone]) continue;
            result[p.sector] = fleet.ids[p.drone];
            droneUsed[p.drone] = 1;
        }
        for (int s = 0; s < S; ++s) {
            if (result[s] >= 0) continue;
            for (int d = 0; d < D; ++d) if (!droneUsed[d]) { result[s] = fleet.ids[d]; droneUsed[d] = 1; break; }
        }
        return result;
    }
}   // namespace ref

// Sector matching now counts shared land over the band cells only, with each previous sector's drone looked up
// once. Over a growing, moving fire, each plan's sectorDrone must equal what the original window scan gives for the
// same plan and fleet. Three kinds of fleet are tried at every step: the previous plan's drones in order; the same
// drones reversed (so fleet index and id differ); and one drone retired plus one the previous plan never saw (so
// some previous sectors have no drone, and one drone has no previous land).
TEST(sector_matching_over_band_matches_window_scan) {
    const std::vector<std::vector<std::uint8_t>> fires = {
        discFire(CX, CY, 40), discFire(CX, CY, 55), bayFire(), discFire(CX + 30, CY - 20, 70), discFire(CX + 30, CY - 20, 90) };
    int compared = 0, differ = 0, matchedByLand = 0;
    MissionPlan prev;
    int nextId = 100;
    for (const auto& fire : fires) {
        DistanceTransform live(GRID_WIDTH, GRID_HEIGHT);
        live.compute(fire);
        Perimeter perim(fire);
        SafePlanner sp(live);
        MissionPlan kept;   // The plan the next fire's fleets follow on from
        for (int variant = 0; variant < 3; ++variant) {
            SweepFleet fleet;
            if (prev.sectorDrone.empty()) {
                for (int i = 0; i < SWEEP_DRONES; ++i) { fleet.ids.push_back(i); fleet.positions.push_back(Vec2(0.5f, 0.5f)); }
            }
            else {
                for (size_t k = 0; k < prev.sectorDrone.size(); ++k) {
                    fleet.ids.push_back(prev.sectorDrone[k]);
                    fleet.positions.push_back(prev.sweepRoutes[k].empty() ? Vec2(0.5f, 0.5f) : prev.sweepRoutes[k].front());
                }
                if (variant == 1) {
                    std::reverse(fleet.ids.begin(), fleet.ids.end());
                    std::reverse(fleet.positions.begin(), fleet.positions.end());
                }
                if (variant == 2 && fleet.ids.size() > 1) {
                    fleet.ids.erase(fleet.ids.begin());             // Retired: its old sector has no drone now
                    fleet.positions.erase(fleet.positions.begin());
                    fleet.ids.push_back(999);                       // New: no land in the previous plan
                    fleet.positions.push_back(Vec2(CX + 0.5f, CY + 0.5f));
                }
                fleet.previous = &prev;
            }
            std::mt19937 rng(2024);
            MissionPlan plan = planMission(live, perim.distance(), sp, rng, fleet);
            const std::vector<int> original = ref::matchSectorsOriginal(plan, plan.sweepRoutes, fleet, plan.bounds);
            differ += plan.sectorDrone != original;
            ++compared;
            // Sectors given a drone already in the fleet (by shared land or distance), so matching really ran
            if (fleet.previous)
                for (int id : plan.sectorDrone) matchedByLand += std::find(fleet.ids.begin(), fleet.ids.end(), id) != fleet.ids.end();
            if (variant == 0) kept = std::move(plan);
        }
        for (int& id : kept.sectorDrone) if (id < 0) id = nextId++;   // As the controller does: new drones get ids
        prev = std::move(kept);
    }
    std::printf("  (%d plans compared, %d differ; %d sectors given a fleet drone)\n", compared, differ, matchedByLand);
    CHECK(differ == 0);
    CHECK(compared == 15);
    CHECK(matchedByLand > 0);
}

// The perimeter-line trace now works out each sample once and slides two row buffers down the window. The
// loops must be exactly those of the original, which worked out all four corners of every square.
TEST(perimeter_trace_matches_per_corner_evaluation) {
    int loopsChecked = 0;
    for (const auto& fire : exactnessFires()) {
        Perimeter perim(fire);
        // The planner's level, the survey standoff, and an odd fractional level
        for (float level : { static_cast<float>(PERIMETER_DISTANCE), static_cast<float>(MIN_FIRE_DISTANCE), 57.3f }) {
            const auto fast = tracePerimeterLine(perim.distance(), level);
            const auto original = ref::tracePerimeterLineOriginal(perim.distance(), level);
            CHECK(fast.size() == original.size());
            if (fast.size() != original.size()) continue;
            for (size_t l = 0; l < fast.size(); ++l) {
                bool same = fast[l].size() == original[l].size();
                for (size_t k = 0; same && k < fast[l].size(); ++k)
                    same = fast[l][k].x == original[l][k].x && fast[l][k].y == original[l][k].y;
                CHECK(same);
                ++loopsChecked;
            }
        }
    }
    std::printf("  (%d loops identical point for point)\n", loopsChecked);
    CHECK(loopsChecked > 0);
}

TEST(sweep_lanes_point_at_a_round_fire) {
    const int r = 30;
    const Vec2 c(400.5f, 400.5f);
    SweepSetup s(discFire(400, 400, r));
    int mainLanes = 0, badDir = 0, badStart = 0, badEnd = 0;
    for (const SweepLane& l : s.plan.lanes) {
        if (l.gapFill) continue;
        ++mainLanes;
        // Direction: from the lane's outer end straight at the fire's centre, to within 3 degrees
        const Vec2 along = l.pts.back() - l.pts.front(), toC = c - l.pts.front();
        const float cosA = (along.x * toC.x + along.y * toC.y)
            / (std::sqrt(along.x * along.x + along.y * along.y) * std::sqrt(toC.x * toC.x + toC.y * toC.y));
        if (cosA < std::cos(3.0f * 3.14159265f / 180.0f)) ++badDir;
        // Runs from the perimeter line (band's outer edge) to the 300 m standoff
        const float rOut = std::sqrt(sqDist(l.pts.front(), c)), rIn = std::sqrt(sqDist(l.pts.back(), c));
        if (std::fabs(rOut - (r + PERIMETER_DISTANCE)) > 2.0f) ++badStart;
        if (std::fabs(rIn - (r + MIN_FIRE_DISTANCE)) > 2.0f) ++badEnd;
    }
    std::printf("  (%d main lanes; off direction %d, off start %d, off end %d)\n", mainLanes, badDir, badStart, badEnd);
    // One lane per 120 m (at most) of perimeter line
    CHECK(mainLanes == static_cast<int>(std::ceil(loopLength(s.plan.perimeterLoops[0]) / LANE_SPACING)));
    CHECK(badDir == 0);
    CHECK(badStart == 0);
    CHECK(badEnd == 0);
}

TEST(sweep_lanes_cover_the_whole_band) {
    SweepSetup round(discFire(400, 400, 30));
    int gapLanes = 0;
    for (const SweepLane& l : round.plan.lanes) gapLanes += l.gapFill;
    CHECK(gapLanes == 0);                      // Convex fire: the main lanes alone see everything
    CHECK(unseenInnerCells(round.plan) == 0);

    SweepSetup bay(bayFire());
    gapLanes = 0;
    for (const SweepLane& l : bay.plan.lanes) gapLanes += l.gapFill;
    std::printf("  (bay: %zu lanes, %d of them gap fill)\n", bay.plan.lanes.size(), gapLanes);
    CHECK(gapLanes > 0);                       // Lanes spread apart in the bay, so gap fill had work to do
    CHECK(unseenInnerCells(bay.plan) == 0);
}

TEST(sweep_sectors_are_exclusive_and_contiguous) {
    SweepSetup s(bayFire());
    const MissionPlan& plan = s.plan;
    const int sectors = static_cast<int>(plan.sweepRoutes.size());

    // Every inner cell has exactly one owner (a valid sector); no cell outside the band has one
    int unowned = 0, strayOwner = 0;
    for (int i = 0; i < N; ++i) {
        if (plan.zone[i] == 1) unowned += plan.owner[i] < 0 || plan.owner[i] >= sectors;
        else strayOwner += plan.owner[i] >= 0;
    }
    CHECK(unowned == 0);
    CHECK(strayOwner == 0);

    // Each sector is one contiguous run of the lane sequence
    std::vector<int> runsOf(sectors, 0);
    for (size_t i = 0; i < plan.lanes.size(); ++i) {
        CHECK(plan.lanes[i].sector >= 0 && plan.lanes[i].sector < sectors);
        if (i == 0 || plan.lanes[i].sector != plan.lanes[i - 1].sector) ++runsOf[plan.lanes[i].sector];
    }
    for (int r : runsOf) CHECK(r == 1);

    // At the perimeter line, the edge lanes of neighbouring sectors are one lane spacing apart (+- 1 cell),
    // i.e. the border is halfway between them and each edge lane is 60 m inside it. Measured on a round fire,
    // where the straight distance between seeds equals the distance along the line (the bay puts a sharp
    // corner in the perimeter line, where the straight distance is shorter).
    SweepSetup round(discFire(400, 400, 30));
    std::vector<const SweepLane*> mainLanes;
    for (const SweepLane& l : round.plan.lanes) if (!l.gapFill) mainLanes.push_back(&l);
    int borders = 0;
    for (size_t i = 0; i < mainLanes.size(); ++i) {
        const SweepLane& a = *mainLanes[i];
        const SweepLane& b = *mainLanes[(i + 1) % mainLanes.size()];
        if (a.sector == b.sector) continue;
        ++borders;
        CHECK(std::fabs(std::sqrt(sqDist(a.seed, b.seed)) - LANE_SPACING) <= 1.0f);
    }
    CHECK(borders == static_cast<int>(round.plan.sweepRoutes.size()));   // One loop: a border per sector
}
