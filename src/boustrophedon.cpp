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
// the few samples that could reach it, nearest first, stopping at the first hit:
//  - A cell whose distance from the run's LINE is more than r + 1 is out of reach of every sample.
//  - Otherwise the samples within reach lie within (r + 1) cells, along the run, of the cell's projection onto
//    it. Samples further along are more than r + 1 away (the margin of 1 cell dwarfs float rounding, which is
//    about 1e-3 cells at these coordinates), so the original test would reject them too.
// The samples that are tested use exactly the original expressions (same float sample point, same distance test),
// so the marked cells are exactly those the disc stamps marked. tests.cpp checks this bit for bit against the
// disc version (lane_stamp_matches_disc_stamps).
void stampLaneCoverage(const SweepLane& lane, const GridRect& area, std::vector<std::uint8_t>& covered) {
    const float r = LANE_COVER_RADIUS, r2 = r * r;
    const float reach = r + 1.0f;   // Search margin: comfortably more than r (see above)
    const int w = area.x1 - area.x0;
    std::vector<Vec2> samples;   // Sample points of the current run (reused between runs)
    for (size_t i = 0; i + 1 < lane.pts.size(); i += 2) {
        const Vec2 a = lane.pts[i], b = lane.pts[i + 1];
        const float len = std::sqrt(sqDist(a, b));
        const int steps = std::max(1, static_cast<int>(std::ceil(len / RAY_STEP)));
        // The run's sample points, worked out once with exactly the original expression
        samples.resize(static_cast<size_t>(steps) + 1);
        for (int s = 0; s <= steps; ++s) samples[s] = a + (b - a) * (static_cast<float>(s) / steps);
        auto sample = [&](int s) { return samples[s]; };

        // A (near) zero-length run has no direction to project onto: stamp its samples as discs
        if (len < 1e-3f) {
            for (int s = 0; s <= steps; ++s) stampDisc(sample(s), area, covered);
            continue;
        }
        const Vec2 dir = (b - a) / len;                    // Unit vector along the run
        const float samplesPerCell = steps / len;          // Sample index advance per cell along the run
        const int span = static_cast<int>(std::ceil(reach * samplesPerCell)) + 1;   // Samples either side to try

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
            for (int x = rx0; x <= rx1; ++x) {
                std::uint8_t& cov = covered[(y - area.y0) * w + (x - area.x0)];
                if (cov) continue;                                     // Already covered (by this or another lane)
                const Vec2 p(x + 0.5f, y + 0.5f), rel = p - a;
                const float along = rel.x * dir.x + rel.y * dir.y;     // Distance along the run from a
                const float across = rel.x * dir.y - rel.y * dir.x;    // Signed distance from the run's line
                if (std::fabs(across) > reach) continue;               // Out of reach of every sample
                if (along < -reach || along > len + reach) continue;   // Beyond either end, out of reach
                // Samples nearest the projection first (most covered cells hit on the first test)
                const int mid = std::clamp(static_cast<int>(std::lround(along * samplesPerCell)), 0, steps);
                for (int k = 0; k <= span; ++k) {
                    const int lo = mid - k, hi = mid + k;
                    if (lo >= 0 && sqDist(p, sample(lo)) <= r2) { cov = 1; break; }
                    if (k > 0 && hi <= steps && sqDist(p, sample(hi)) <= r2) { cov = 1; break; }
                    if (lo < 0 && hi > steps) break;                   // No samples left on either side
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
// Ownership Map
// ----------------------------------------------------------------------------
std::vector<int> buildOwnerMap(const std::vector<SweepLane>& lanes, const std::vector<std::uint8_t>& zone,
                               const GridRect& area) {
    std::vector<int> owner(static_cast<size_t>(GRID_WIDTH) * GRID_HEIGHT, -1);
    buildOwnerMap(lanes, zone, area, owner);
    return owner;
}

// Fills a caller-provided whole-grid array (all -1 on entry); only cells inside `area` are written
void buildOwnerMap(const std::vector<SweepLane>& lanes, const std::vector<std::uint8_t>& zone,
                   const GridRect& area, std::vector<int>& owner) {
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
                if (zone[idx] != 1 || owner[idx] >= 0) continue;
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
    const std::uint8_t* const zp = zone.data();
    int* const op = owner.data();
    for (size_t head = 0; head < queue.size(); ++head) {
        const int idx = queue[head];
        const int x = idx % GRID_WIDTH, y = idx / GRID_WIDTH;
        const int sector = op[idx];
        if (x > area.x0 && x < area.x1 - 1 && y > area.y0 && y < area.y1 - 1) {
            // Interior cell: all eight neighbours are inside the area
            for (int k = 0; k < NUM_DIRS; ++k) {
                const int ni = idx + offset[k];
                if (zp[ni] != 1 || op[ni] >= 0) continue;
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
            if (zp[ni] != 1 || op[ni] >= 0) continue;
            op[ni] = sector;
            queue.push_back(ni);
        }
    }

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
