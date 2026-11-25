#include "boustrophedon.h"
#include <algorithm>
#include <cmath>
#include <limits>

// ============================================================================
// Sweep Lanes at Right Angles to the Fire
//
// Replaces the old north-south lanes per Voronoi region. Lanes now start on the perimeter line and run
// straight in towards the fire. For an exact Euclidean distance field the lines of steepest descent are
// straight lines to the nearest fire point, so one direction sample per seed gives a straight lane at right
// angles to the nearest fire edge.
// ============================================================================

namespace {
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

    // Disc stamp at one sample point c: marks every cell of `area` whose centre is within LANE_COVER_RADIUS of c.
    // This is the original definition of coverage; the per-cell version in stampLane must agree with it exactly.
    void stampDisc(Vec2 c, const GridRect& area, std::vector<std::uint8_t>& covered) {
        const float r = LANE_COVER_RADIUS, r2 = r * r;
        const int w = area.x1 - area.x0;
        const int x0 = std::max(area.x0, static_cast<int>(std::floor(c.x - r)));
        const int x1 = std::min(area.x1 - 1, static_cast<int>(std::floor(c.x + r)));
        const int y0 = std::max(area.y0, static_cast<int>(std::floor(c.y - r)));
        const int y1 = std::min(area.y1 - 1, static_cast<int>(std::floor(c.y + r)));
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x)
                if (sqDist(Vec2(x + 0.5f, y + 0.5f), c) <= r2) covered[(y - area.y0) * w + (x - area.x0)] = 1;
    }
}

// ----------------------------------------------------------------------------
// Lane coverage stamp
// ----------------------------------------------------------------------------
// Marks every cell (inside area, as a local bitmap) whose centre is within LANE_COVER_RADIUS of one of the
// lane's samples: each run a -> b is sampled every RAY_STEP (sample s of `steps` is a + (b - a) * (s / steps)),
// and a cell is covered when it is within the radius of a sample. A cell within the radius of a sample is within
// the radius of the run, so a marked cell is always really covered.
//
// Coverage used to be stamped as a disc (a 14 x 14 box test) at every sample. Neighbouring samples are only half
// a cell apart, so each cell was tested against ~28 samples. Now each cell near the run is tested only against
// the few samples that could reach it, stopping at the first hit:
//  - Every sample lies on the run's line (to within float rounding, about 1e-4 cells at these coordinates). So a
//    cell whose distance `across` from the line is more than r + MARGIN is out of reach of every sample.
//    MARGIN (0.1 cell) is about a thousand times the rounding, so the original float test would reject those
//    samples too.
//  - A cell level with the run and within sqrt((r - MARGIN)^2 - (RAY_STEP / 2)^2) of its line is certainly
//    covered: the nearest sample is at most half a step along, so it is within r - MARGIN of the cell and passes
//    the original test. Most covered cells are marked this way, with no sample test at all.
//  - Otherwise a sample can only be within r of the cell if it lies within h = sqrt((r + MARGIN)^2 - across^2)
//    of the cell's projection onto the run (Pythagoras; MARGIN again absorbs rounding). Sample s sits s / (steps
//    per cell) along the run, so those samples form one index range, widened by one sample each side. Only that
//    range is tested. (Until 2026-11-23 every cell within r + 1 of the line tried all samples within r + 1 of its
//    projection, about 30, so the cells just outside the swath, which no sample reaches, cost the most.)
// The samples that are tested use exactly the original expressions (same float sample point, same distance test),
// and the samples skipped provably fail that test, so the marked cells are exactly those the disc stamps marked.
// tests.cpp checks this bit for bit against the disc version (lane_stamp_matches_disc_stamps, and
// lane_stamp_matches_disc_stamps_on_the_edge for cells right on the swath's edge).
void stampLaneCoverage(const SweepLane& lane, const GridRect& area, std::vector<std::uint8_t>& covered) {
    const float r = LANE_COVER_RADIUS, r2 = r * r;
    constexpr float MARGIN = 0.1f;          // Allowance for float rounding (see above), far larger than it
    const float reach = r + MARGIN;         // No sample is within r of a cell further than this from the line
    const float reach2 = reach * reach;
    // Cells within sqrt(inner2) of the line (and level with the run) are certainly covered (see the loop below)
    const float inner2 = (r - MARGIN) * (r - MARGIN) - 0.25f * RAY_STEP * RAY_STEP;
    const int w = area.x1 - area.x0;
    std::vector<Vec2> samples;   // Sample points of the current run (reused between runs)
    for (size_t i = 0; i + 1 < lane.pts.size(); i += 2) {
        const Vec2 a = lane.pts[i], b = lane.pts[i + 1];
        const float len = std::sqrt(sqDist(a, b));
        const int steps = std::max(1, static_cast<int>(std::ceil(len / RAY_STEP)));
        // The run's sample points, worked out once with exactly the original expression
        samples.resize(static_cast<size_t>(steps) + 1);
        for (int s = 0; s <= steps; ++s) samples[s] = a + (b - a) * (static_cast<float>(s) / steps);

        // A (near) zero-length run has no direction to project onto: stamp its samples as discs
        if (len < 1e-3f) {
            for (int s = 0; s <= steps; ++s) stampDisc(samples[s], area, covered);
            continue;
        }
        const Vec2 dir = (b - a) / len;                    // Unit vector along the run
        const float samplesPerCell = steps / len;          // Sample index advance per cell along the run

        // Cells whose centres could be within reach: the run's bounding box grown by `reach`, clipped to area
        const int x0 = std::max(area.x0, static_cast<int>(std::floor(std::min(a.x, b.x) - reach)));
        const int x1 = std::min(area.x1 - 1, static_cast<int>(std::floor(std::max(a.x, b.x) + reach)));
        const int y0 = std::max(area.y0, static_cast<int>(std::floor(std::min(a.y, b.y) - reach)));
        const int y1 = std::min(area.y1 - 1, static_cast<int>(std::floor(std::max(a.y, b.y) + reach)));
        for (int y = y0; y <= y1; ++y) {
            // On this row, only cells whose distance from the run's line could be within `reach` need looking at.
            // That distance is linear in x, so they form one interval; work it out (one cell wider each side, so
            // rounding here can't drop a cell: the exact per-cell test below still decides).
            int rx0 = x0, rx1 = x1;
            if (std::fabs(dir.y) > 1e-6f) {
                const float base = (y + 0.5f - a.y) * dir.x;          // across = (x + 0.5 - a.x) * dir.y - base
                const float xa = a.x - 0.5f + (base - reach) / dir.y, xb = a.x - 0.5f + (base + reach) / dir.y;
                rx0 = std::max(x0, static_cast<int>(std::floor(std::min(xa, xb))) - 1);
                rx1 = std::min(x1, static_cast<int>(std::ceil(std::max(xa, xb))) + 1);
            }
            std::uint8_t* const row = covered.data() + static_cast<size_t>(y - area.y0) * w - area.x0;
            for (int x = rx0; x <= rx1; ++x) {
                std::uint8_t& cov = row[x];
                if (cov) continue;                                     // Already covered (by this or another lane)
                const Vec2 p(x + 0.5f, y + 0.5f), rel = p - a;
                const float across = rel.x * dir.y - rel.y * dir.x;    // Signed distance from the run's line
                const float acrossSq = across * across;
                if (acrossSq > reach2) continue;                       // Out of reach of every sample
                const float along = rel.x * dir.x + rel.y * dir.y;     // Distance along the run from a
                // Well inside the swath: the nearest sample is at most half a step (RAY_STEP / 2) along from the
                // cell's projection, so within sqrt(across^2 + (RAY_STEP / 2)^2) <= r - MARGIN of the cell, and
                // the original float test certainly passes for it. Mark the cell without testing any sample.
                if (acrossSq <= inner2 && along >= 0.0f && along <= len) { cov = 1; continue; }
                // Near the swath's edge or past an end: the samples that could be within r are those within h
                // of the projection (see the header), as an index range widened by one sample each side and
                // clipped to the run's samples
                const float h = std::sqrt(reach2 - acrossSq);
                const float sLo = (along - h) * samplesPerCell, sHi = (along + h) * samplesPerCell;
                if (sHi < -1.0f || sLo > steps + 1.0f) continue;       // Beyond either end, out of reach
                const int lo = std::max(0, static_cast<int>(std::floor(sLo)) - 1);
                const int hi = std::min(steps, static_cast<int>(std::ceil(sHi)) + 1);
                // Try the sample nearest the projection first, then outwards (the order only affects how soon
                // a hit stops the search, not whether there is one)
                const int mid = std::clamp(static_cast<int>(std::floor(along * samplesPerCell + 0.5f)), lo, hi);
                for (int k = 0; mid - k >= lo || mid + k <= hi; ++k) {
                    if (mid - k >= lo && sqDist(p, samples[mid - k]) <= r2) { cov = 1; break; }
                    if (k > 0 && mid + k <= hi && sqDist(p, samples[mid + k]) <= r2) { cov = 1; break; }
                }
            }
        }
    }
}

namespace {
}

// ----------------------------------------------------------------------------
// Lay Lanes Along the Perimeter Line
// ----------------------------------------------------------------------------
// For each loop: n = ceil(length / LANE_SPACING) seeds, evenly spaced (so never more than 120 m apart and
// the spacing has no short leftover at the end of the loop). Sector borders fall halfway between two lanes,
// so each sector's first lane is half a spacing (60 m) inside its border.
std::vector<SweepLane> layLanes(const std::vector<std::vector<Vec2>>& loops, const DistanceTransform& perimEdt,
                                const std::vector<std::uint8_t>& zone) {
    std::vector<SweepLane> lanes;
    for (int li = 0; li < static_cast<int>(loops.size()); ++li) {
        const std::vector<Vec2>& loop = loops[li];
        const size_t m = loop.size();
        float total = 0.0f;
        for (size_t i = 0; i < m; ++i) total += std::sqrt(sqDist(loop[i], loop[(i + 1) % m]));
        if (total <= 0.0f) continue;
        const int n = std::max(1, static_cast<int>(std::ceil(total / LANE_SPACING)));
        const float spacing = total / n;

        // Walk the loop once, dropping a seed every `spacing` cells of arc length
        size_t seg = 0;
        float segStart = 0.0f;   // Arc length at the start of segment `seg`
        for (int k = 0; k < n; ++k) {
            const float want = k * spacing;
            float segLen = std::sqrt(sqDist(loop[seg], loop[(seg + 1) % m]));
            while (segStart + segLen < want && seg + 1 < m) {
                segStart += segLen;
                ++seg;
                segLen = std::sqrt(sqDist(loop[seg], loop[(seg + 1) % m]));
            }
            const Vec2 a = loop[seg], b = loop[(seg + 1) % m];
            const float u = segLen > 1e-6f ? (want - segStart) / segLen : 0.0f;
            const Vec2 seed = a + (b - a) * u;

            // Loops run anticlockwise on screen with the fire on the left, so in grid coordinates (y down) the
            // inward normal of travel direction (tx, ty) is (ty, -tx). Used only where the gradient fails.
            const Vec2 tangent = b - a;
            const Vec2 dir = downhill(perimEdt, seed, Vec2(tangent.y, -tangent.x));

            // Start a little outside the line so the lane reaches the band's outer edge too, and walk in at most
            // twice the band's width (the standoff is normally reached after PERIMETER - MIN_FIRE cells)
            SweepLane lane;
            lane.pts = bandRuns(perimEdt, zone, seed, dir, -2.0f, 2.0f * PERIMETER_DISTANCE, true);
            if (lane.pts.empty()) continue;   // Fire covers the band here
            lane.seed = seed;
            lane.loop = li;
            lanes.push_back(std::move(lane));
        }
    }
    return lanes;
}

// ----------------------------------------------------------------------------
// Gap Fill in Concave Bays
// ----------------------------------------------------------------------------
int fillCoverageGaps(std::vector<SweepLane>& lanes, const DistanceTransform& perimEdt,
                     const std::vector<std::uint8_t>& zone, const GridRect& area) {
    if (area.empty() || lanes.empty()) return 0;
    const int w = area.x1 - area.x0, h = area.y1 - area.y0;
    std::vector<std::uint8_t> covered(static_cast<size_t>(w) * h, 0);
    for (const SweepLane& l : lanes) stampLaneCoverage(l, area, covered);

    constexpr int MAX_PASSES = 6;
    int added = 0;

    // Cells that may start a patch, as local indices in increasing order. The first pass needs every unseen
    // inner cell in the window, found with one row-by-row scan. A later pass only needs the cells of the
    // previous pass's patches: lanes only ever ADD coverage, so a cell unseen now was unseen in the previous
    // pass too, and was therefore in one of its patches. Looking at those cells in the same increasing order
    // starts every flood fill from the same cell as a scan of the whole window, and the flood fill itself is
    // unchanged, so the patches (and the order of the cells in each) are exactly the same as before.
    std::vector<int> candidates;
    for (int y = 0; y < h; ++y) {
        const std::uint8_t* zoneRow = zone.data() + static_cast<size_t>(area.y0 + y) * GRID_WIDTH + area.x0;
        for (int x = 0; x < w; ++x) {
            if (zoneRow[x] == 1 && !covered[y * w + x]) candidates.push_back(y * w + x);
        }
    }
    // Flood-fill marks, one byte per window cell: equal to pass + 1 once the cell is in a patch of this pass.
    // Allocated once; a new pass uses a new mark value instead of clearing the array (MAX_PASSES < 255).
    std::vector<std::uint8_t> inPatch(static_cast<size_t>(w) * h, 0);

    for (int pass = 0; pass < MAX_PASSES; ++pass) {
        const std::uint8_t mark = static_cast<std::uint8_t>(pass + 1);
        // --- Find patches of unseen inner cells (8-connected flood fill) ---
        std::vector<std::vector<int>> patches;   // Local indices of the cells in each patch
        for (int i : candidates) {
            // Candidates are inner cells; skip the ones a new lane now covers or an earlier patch took
            if (covered[i] || inPatch[i] == mark) continue;
            const int id = static_cast<int>(patches.size());
            patches.emplace_back();
            std::vector<int> stack{ i };
            inPatch[i] = mark;
            while (!stack.empty()) {
                const int c = stack.back();
                stack.pop_back();
                patches[id].push_back(c);
                const int cx = c % w, cy = c / w;
                for (int k = 0; k < NUM_DIRS; ++k) {
                    const int nx = cx + DX[k], ny = cy + DY[k];
                    if (nx < 0 || nx >= w || ny < 0 || ny >= h) continue;
                    const int ni = ny * w + nx;
                    if (covered[ni] || inPatch[ni] == mark) continue;
                    if (zone[(area.y0 + ny) * GRID_WIDTH + area.x0 + nx] != 1) continue;
                    inPatch[ni] = mark;
                    stack.push_back(ni);
                }
            }
        }
        if (patches.empty()) break;   // Everything is seen

        // Next pass's candidates: this pass's patch cells, back in increasing order
        candidates.clear();
        for (const std::vector<int>& patch : patches) candidates.insert(candidates.end(), patch.begin(), patch.end());
        std::sort(candidates.begin(), candidates.end());

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
            stampLaneCoverage(g, area, covered);
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

// ----------------------------------------------------------------------------
// Sector Route and its Estimated Lap
// ----------------------------------------------------------------------------
std::vector<Vec2> sectorWaypoints(const std::vector<SweepLane>& lanes, int first, int last) {
    std::vector<Vec2> wp;
    for (int i = first; i <= last; ++i) {
        const std::vector<Vec2>& p = lanes[i].pts;
        const bool inward = (i - first) % 2 == 0;   // Alternate: in, out, in, ...
        auto push = [&](Vec2 v) { if (wp.empty() || sqDist(wp.back(), v) > 1e-6f) wp.push_back(v); };
        if (inward) for (const Vec2& v : p) push(v);
        else for (auto it = p.rbegin(); it != p.rend(); ++it) push(*it);
    }
    return wp;
}

float estimateSectorLength(const std::vector<SweepLane>& lanes, int first, int last) {
    // One builder walk over the lanes; SectorLapBuilder holds the arithmetic so the two can never disagree
    SectorLapBuilder b(lanes, first);
    while (b.last() < last) b.addNext();
    return b.lap();
}

// ----------------------------------------------------------------------------
// Incremental Sector Lap (SectorLapBuilder)
// ----------------------------------------------------------------------------
// The lap is: first lane's length, then for each further lane the straight step from the previous lane's end
// to its entry plus its length, then the straight return leg to the first lane's entry. Lanes alternate
// direction (inward, outward, ...) counted from the sector's first lane. Keeping the running open length means
// adding a lane, or asking what the lap would be with one more, costs one lane's worth of work.
SectorLapBuilder::SectorLapBuilder(const std::vector<SweepLane>& lanes, int first)
    : lanes_(lanes), first_(first), last_(first) {
    const std::vector<Vec2>& p = lanes[first].pts;
    start_ = p.front();                      // The first lane is flown inward
    end_ = p.back();
    len_ += lanes[first].length();           // 0 + length, as estimateSectorLength always did
}

void SectorLapBuilder::extended(int i, float& len, Vec2& end) const {
    const std::vector<Vec2>& p = lanes_[i].pts;
    const bool inward = (i - first_) % 2 == 0;
    const Vec2 a = inward ? p.front() : p.back(), b = inward ? p.back() : p.front();
    len += std::sqrt(sqDist(end, a));        // Step from the previous lane (kept as two separate additions,
    len += lanes_[i].length();               // in this order, so results match the old loop bit for bit)
    end = b;
}

float SectorLapBuilder::lap() const { return len_ + std::sqrt(sqDist(end_, start_)); }

float SectorLapBuilder::lapWithNext() const {
    float len = len_;
    Vec2 end = end_;
    extended(last_ + 1, len, end);
    return len + std::sqrt(sqDist(end, start_));   // Return leg to close the loop
}

void SectorLapBuilder::addNext() { extended(++last_, len_, end_); }

// ----------------------------------------------------------------------------
// Memoised Sector Laps (SectorLapTable)
// ----------------------------------------------------------------------------
// The sector split binary-searches a lap limit, and every step of that search re-runs a greedy fill that grows
// sectors lane by lane from the same few first lanes. Laps are pure functions of (first, last) and the lanes, so
// each is worked out once here and then only looked up. A row per first lane holds one SectorLapBuilder (the
// only place the arithmetic lives) and the laps it has produced, so every value is bit for bit the one the
// builder, and so estimateSectorLength, gives.
SectorLapTable::SectorLapTable(const std::vector<SweepLane>& lanes) : lanes_(lanes) { reset(); }

void SectorLapTable::reset() {
    rows_.clear();
    rows_.resize(lanes_.size());   // All empty: rows are created on first use
}

SectorLapTable::Row& SectorLapTable::row(int first) {
    std::unique_ptr<Row>& r = rows_[static_cast<size_t>(first)];
    if (!r) {
        r = std::make_unique<Row>(lanes_, first);
        r->laps.push_back(r->builder.lap());                          // laps[0]: the first lane alone
        r->runMax.push_back(-std::numeric_limits<float>::infinity()); // Unused (the first lane is always taken)
    }
    return *r;
}

void SectorLapTable::grow(Row& r) {
    r.builder.addNext();
    const float l = r.builder.lap();   // Same value lapWithNext() gave before addNext() (tested bit for bit)
    r.laps.push_back(l);
    // A NaN lap (only possible with non-finite lane coordinates) fails every `<= limit` test, so for the
    // running maximum it counts as +infinity, which fails every finite limit too
    const float key = std::isnan(l) ? std::numeric_limits<float>::infinity() : l;
    r.runMax.push_back(std::max(r.runMax.back(), key));
}

float SectorLapTable::lap(int first, int last) {
    Row& r = row(first);
    const size_t k = static_cast<size_t>(last - first);
    while (r.laps.size() <= k) grow(r);
    return r.laps[k];
}

// The original greedy fill was:
//     SectorLapBuilder b(lanes, first);
//     while (b.last() + 1 <= bound && b.lapWithNext() <= limit) b.addNext();
//     return b.last();
// i.e. it stops just before the first lane j in (first, bound] whose lap(first, j) is not <= limit. With a finite
// limit, "some lap up to j fails" is exactly "the running maximum at j is > limit", and the running maximum never
// falls, so the stopping lane among the laps already worked out is found by binary search. Past those, the row
// grows one lane at a time, as the original did, until a lap fails or the bound is reached.
int SectorLapTable::greedyLast(int first, int bound, float limit) {
    if (bound <= first) return first;
    Row& r = row(first);
    const int need = bound - first;   // Highest k = last - first the fill may reach

    if (std::isfinite(limit)) {
        const int have = std::min(need, static_cast<int>(r.laps.size()) - 1);   // laps[1..have] known
        const auto b = r.runMax.begin() + 1, e = r.runMax.begin() + 1 + have;
        const auto it = std::upper_bound(b, e, limit);                         // First running max > limit
        if (it != e) return first + static_cast<int>(it - r.runMax.begin()) - 1;
        for (int k = have + 1; k <= need; ++k) {   // have < need here, so laps.size() == k before grow
            grow(r);
            if (!(r.laps[static_cast<size_t>(k)] <= limit)) return first + k - 1;
        }
        return bound;
    }

    // Infinite or NaN limit (never produced by the planner): the plain lane-by-lane test, for exactness
    for (int k = 1; k <= need; ++k) {
        if (r.laps.size() <= static_cast<size_t>(k)) grow(r);
        if (!(r.laps[static_cast<size_t>(k)] <= limit)) return first + k - 1;
    }
    return bound;
}

// ----------------------------------------------------------------------------
// Ownership Map
// ----------------------------------------------------------------------------
std::vector<int> buildOwnerMap(const std::vector<SweepLane>& lanes, const std::vector<std::uint8_t>& zone,
                               const GridRect& area, long long bandCells) {
    std::vector<int> owner(static_cast<size_t>(GRID_WIDTH) * GRID_HEIGHT, -1);
    buildOwnerMap(lanes, zone, area, owner, bandCells);
    return owner;
}

// Fills a caller-provided whole-grid array (all -1 on entry); only cells inside `area` are written. The zone map
// is read-only here, so the flood runs on a copy of it (the planner, which owns its zone map, calls
// buildOwnerMapUsingZone directly and saves the copy).
void buildOwnerMap(const std::vector<SweepLane>& lanes, const std::vector<std::uint8_t>& zone,
                   const GridRect& area, std::vector<int>& owner, long long bandCells) {
    std::vector<std::uint8_t> scratch = zone;
    buildOwnerMapUsingZone(lanes, scratch, area, owner, bandCells);
}

// ----------------------------------------------------------------------------
// Ownership Map, flooding through the caller's zone map
// ----------------------------------------------------------------------------
// The same map as buildOwnerMap, built the same way, but the flood keeps track of which band cells are still
// free in the zone map itself instead of reading the owner array to find out.
//
// Why: the flood tries 8 neighbours of every band cell, and each try asked two questions, "is it a band cell?"
// (one byte of `zone`) and "has it an owner yet?" (four bytes of `owner`, a second array). Marking each band
// cell in `zone` as it is claimed (ZONE_CLAIMED) turns both into one question about one byte: a neighbour is
// free exactly when its zone is still 1. That halves the arrays the flood reads, and the 4-byte owner array is
// now only written, once per cell.
//
// Exactness: a cell is marked at the very moment it is given an owner (seeding or flood), so "zone == 1" here
// means precisely "zone == 1 and owner < 0" in the original test. The cells, their order in the queue and the
// sector each gets are therefore unchanged. Every marked cell is in the queue, so walking the queue afterwards
// puts every mark back to 1, and the zone map is exactly as it was on entry before the final scan reads it.
void buildOwnerMapUsingZone(const std::vector<SweepLane>& lanes, std::vector<std::uint8_t>& zone,
                            const GridRect& area, std::vector<int>& owner, long long bandCells) {
    constexpr std::uint8_t ZONE_CLAIMED = 3;   // Band cell already given an owner (only while this function runs)
    std::vector<int> queue;

    // Seed the search with every cell a lane run passes over
    for (const SweepLane& l : lanes) {
        for (size_t i = 0; i + 1 < l.pts.size(); i += 2) {
            const Vec2 a = l.pts[i], b = l.pts[i + 1];
            const int steps = std::max(1, static_cast<int>(std::ceil(std::sqrt(sqDist(a, b)) / RAY_STEP)));
            for (int s = 0; s <= steps; ++s) {
                const Vec2 p = a + (b - a) * (static_cast<float>(s) / steps);
                const int x = static_cast<int>(std::floor(p.x)), y = static_cast<int>(std::floor(p.y));
                if (x < area.x0 || x >= area.x1 || y < area.y0 || y >= area.y1) continue;
                const int idx = y * GRID_WIDTH + x;
                if (zone[idx] != 1) continue;   // Not a band cell, or already claimed
                zone[idx] = ZONE_CLAIMED;
                owner[idx] = l.sector;
                queue.push_back(idx);
            }
        }
    }

    // Breadth-first through the band (8-connected): each cell takes the sector of the lane that reaches it first.
    // Cells are taken off the queue in order and their neighbours tried in the fixed DX/DY order, so which lane
    // reaches a cell first never depends on anything else. Two exact speed-ups:
    //  - Cells at least one cell inside the area (nearly all of them) need no bounds test on their neighbours,
    //    so those cells use a loop over precomputed index offsets; only cells on the area's edge test bounds.
    //  - Every band cell enters the queue at most once, so the queue never needs more room than the area has
    //    cells: reserving that up front saves the repeated reallocate-and-copy as it grows.
    queue.reserve(static_cast<size_t>(std::max(0, area.x1 - area.x0)) * std::max(0, area.y1 - area.y0));
    int offset[NUM_DIRS];
    for (int k = 0; k < NUM_DIRS; ++k) offset[k] = DY[k] * GRID_WIDTH + DX[k];
    std::uint8_t* const zp = zone.data();
    int* const op = owner.data();
    for (size_t head = 0; head < queue.size(); ++head) {
        const int idx = queue[head];
        const int x = idx % GRID_WIDTH, y = idx / GRID_WIDTH;
        const int sector = op[idx];
        if (x > area.x0 && x < area.x1 - 1 && y > area.y0 && y < area.y1 - 1) {
            // Interior cell: all eight neighbours are inside the area
            for (int k = 0; k < NUM_DIRS; ++k) {
                const int ni = idx + offset[k];
                if (zp[ni] != 1) continue;   // Not a band cell, or already claimed
                zp[ni] = ZONE_CLAIMED;
                op[ni] = sector;
                queue.push_back(ni);
            }
            continue;
        }
        // Edge cell: the original bounds-checked test, same neighbour order
        for (int k = 0; k < NUM_DIRS; ++k) {
            const int nx = x + DX[k], ny = y + DY[k];
            if (nx < area.x0 || nx >= area.x1 || ny < area.y0 || ny >= area.y1) continue;
            const int ni = ny * GRID_WIDTH + nx;
            if (zp[ni] != 1) continue;   // Not a band cell, or already claimed
            zp[ni] = ZONE_CLAIMED;
            op[ni] = sector;
            queue.push_back(ni);
        }
    }

    // Put the zone map back as it was: every claimed cell is in the queue, and was a band cell (zone 1)
    for (const int idx : queue) zp[idx] = 1;

    // Every band cell enters the queue exactly once, when it is given an owner. If the queue holds as many cells
    // as the area has band cells (the caller's count), none is left without an owner and the scan below would
    // change nothing, so it is skipped. That is the usual case: cut-off slivers are rare.
    if (bandCells >= 0 && static_cast<long long>(queue.size()) == bandCells) return;

    // Band cells cut off from every lane (a sliver beyond a gap in the band): nearest lane in a straight line
    for (int y = area.y0; y < area.y1; ++y) {
        for (int x = area.x0; x < area.x1; ++x) {
            const int idx = y * GRID_WIDTH + x;
            if (zone[idx] != 1 || owner[idx] >= 0 || lanes.empty()) continue;
            const Vec2 p(x + 0.5f, y + 0.5f);
            float bd = std::numeric_limits<float>::max();
            for (const SweepLane& l : lanes) {
                const float d = laneDist(l, p);
                if (d < bd) { bd = d; owner[idx] = l.sector; }
            }
        }
    }
}
