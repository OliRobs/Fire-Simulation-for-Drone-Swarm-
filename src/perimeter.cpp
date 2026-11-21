#include "perimeter.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include "viewport.h"   // maskBounds

// ============================================================================
// Fire Perimeter Tracking & Containment System
//
// Manages the "committed" containment perimeter around the expanding wildfire.
// Instead of advancing continuously with every minor fire spark, the perimeter 
// updates discretely in local "sections" when the live fire encroaches within 
// a specific trigger threshold distance.
// ============================================================================

// Constructor: Initializes committed fire mask and computes the baseline 
// Euclidean Distance Transform (EDT) around the initial fire shape.
Perimeter::Perimeter(const std::vector<std::uint8_t>& initialFire)
    : committed_(initialFire), trigger_(initialFire.size(), 0), lastLive_(initialFire),
    // Range limits (distance_transform.h): the committed-fire distance is read out to SEARCH_DISTANCE + 0.5
    // (mission zones; the ring and the trigger distance are nearer), the trigger distance only out to
    // PERIMETER_SECTION_RADIUS, so neither transform needs the whole grid.
    perimEdt_(GRID_WIDTH, GRID_HEIGHT, SEARCH_DISTANCE + 1),
    triggerEdt_(GRID_WIDTH, GRID_HEIGHT, PERIMETER_SECTION_RADIUS + 1) {
    // Initial distance field computation around committed fire boundaries
    perimEdt_.compute(committed_);
    // Record where the starting fire is, so later updates know where committed_ can be non-zero
    int x0, y0, x1, y1;
    if (maskBounds(committed_, GRID_WIDTH, GRID_HEIGHT, x0, y0, x1, y1)) fireSeen_ = GridRect{ x0, y0, x1 + 1, y1 + 1 };
}

// ----------------------------------------------------------------------------
// Update without a known fire rectangle: find it, then do the normal update
// ----------------------------------------------------------------------------
bool Perimeter::update(const std::vector<std::uint8_t>& liveFire) {
    int x0, y0, x1, y1;
    if (!maskBounds(liveFire, GRID_WIDTH, GRID_HEIGHT, x0, y0, x1, y1)) return update(liveFire, GridRect{});
    return update(liveFire, GridRect{ x0, y0, x1 + 1, y1 + 1 });
}

// ----------------------------------------------------------------------------
// Update Containment Perimeter
// ----------------------------------------------------------------------------
// Evaluates live fire expansion against the committed perimeter distance field.
// Returns `true` if any section of the committed perimeter was pushed forward.
//
// Only cells that have newly caught fire since the last call can trigger an update, so only those are tested.
// The reason: when a call ends, no live cell is within the trigger distance any more.
//  - If the call returned false, no live cell triggered (that is what false means).
//  - If it returned true, every triggering cell is within PERIMETER_SECTION_RADIUS of itself, so it was
//    committed (distance 0 afterwards); every other live cell only got nearer to the committed fire.
// The committed field only changes inside update(), so at the next call a cell that was live last time still
// does not trigger, and the cells to test are those live now that were not live then. They are found by
// comparing the mask with a copy of the last one (lastLive_), row by row (memcmp), which is much cheaper than
// testing every cell of the fire's area against the distance field on every frame. The triggering set, and so
// everything computed from it, is exactly the one the full test would find. This holds whether or not the
// mask only grows: a cell that went out and caught fire again counts as new and is tested again, which is
// harmless (the committed fire only grows, so a cell that was once too near to trigger stays too near).
bool Perimeter::update(const std::vector<std::uint8_t>& liveFire, const GridRect& fireBounds) {
    // Trigger Threshold: A live fire cell triggers a perimeter update once it 
    // encroaches within PERIMETER_UPDATE_DISTANCE of the current perimeter ring.
    constexpr double triggerSq = static_cast<double>(PERIMETER_DISTANCE - PERIMETER_UPDATE_DISTANCE)
        * (PERIMETER_DISTANCE - PERIMETER_UPDATE_DISTANCE);

    // Section Expansion Radius: Maximum influence distance around triggering cells
    constexpr double sectionSq = static_cast<double>(PERIMETER_SECTION_RADIUS) * PERIMETER_SECTION_RADIUS;

    // Every cell that is (or ever was) live fire, or committed, lies inside 'area', so the loops below only
    // visit it: outside it liveFire, lastLive_ and committed_ are all 0. Using the union with the earlier
    // rectangles (not just this one) also catches cells that were live last time but have gone out since.
    fireSeen_ = fireSeen_.unite(fireBounds);
    const GridRect area = fireSeen_;

    // --- Phase 1: Find the newly lit cells that breach the trigger distance ---
    newTriggers_.clear();
    const std::size_t rowBytes = static_cast<std::size_t>(area.x1 - area.x0);
    for (int y = area.y0; y < area.y1 && rowBytes > 0; ++y) {
        const int rowStart = y * GRID_WIDTH + area.x0;
        const std::uint8_t* now = liveFire.data() + rowStart;
        std::uint8_t* last = lastLive_.data() + rowStart;
        if (std::memcmp(now, last, rowBytes) == 0) continue;   // Row unchanged: nothing new can trigger here
        for (std::size_t k = 0; k < rowBytes; ++k) {
            // Newly lit (live now, not live last call) and still at least the trigger distance from the
            // committed fire: exactly the full test, limited to the cells that can pass it
            if (now[k] && !last[k] && perimEdt_.squaredDistance(rowStart + static_cast<int>(k)) >= triggerSq) {
                newTriggers_.push_back(rowStart + static_cast<int>(k));
            }
        }
        std::memcpy(last, now, rowBytes);   // Remember this row for the next call
    }

    // No cells breached the trigger distance; perimeter remains unchanged
    if (newTriggers_.empty()) return false;

    // --- Phase 2: Pull the committed perimeter forward around triggered sections ---
    // Compute distance field from the triggering cells (all inside 'area', so no grid scan is needed).
    // trigger_ is all 0 between calls, so after marking these cells it holds exactly the triggering set.
    for (int idx : newTriggers_) trigger_[idx] = 1;
    triggerEdt_.compute(trigger_, area);
    for (int idx : newTriggers_) trigger_[idx] = 0;   // Back to all 0 for the next call

    for (int y = area.y0; y < area.y1; ++y) {
        for (int idx = y * GRID_WIDTH + area.x0, end = y * GRID_WIDTH + area.x1; idx < end; ++idx) {
            // Expand committed fire mask to include live fire cells within the section radius
            if (liveFire[idx] && !committed_[idx] && triggerEdt_.squaredDistance(idx) <= sectionSq) {
                committed_[idx] = 1;
            }
        }
    }

    // Recompute committed distance field following perimeter updates (committed_ lies inside 'area')
    perimEdt_.compute(committed_, area);
    return true;
}

// ----------------------------------------------------------------------------
// Extract Perimeter Ring Cells
// ----------------------------------------------------------------------------
// Returns a list of grid cell indices forming the visual/functional perimeter 
// ring band located at standard `PERIMETER_DISTANCE` away from committed fire.
std::vector<int> Perimeter::ringCells() const {
    constexpr double ringLo = (PERIMETER_DISTANCE - 0.5) * (PERIMETER_DISTANCE - 0.5);
    constexpr double ringHi = (PERIMETER_DISTANCE + 0.5) * (PERIMETER_DISTANCE + 0.5);

    // The ring is nearer than the transform's range (SEARCH_DISTANCE), so it lies inside the transform's
    // window: only that needs scanning, in the same row-major order as a whole-grid scan
    const GridRect win = perimEdt_.window();
    std::vector<int> ring;
    for (int y = win.y0; y < win.y1; ++y) {
        for (int idx = y * GRID_WIDTH + win.x0, end = y * GRID_WIDTH + win.x1; idx < end; ++idx) {
            const double d2 = perimEdt_.squaredDistance(idx);
            // Collect cells residing within the thin ring shell radius [PERIMETER_DISTANCE ± 0.5]
            if (d2 >= ringLo && d2 < ringHi) {
                ring.push_back(idx);
            }
        }
    }
    return ring;
}

// ============================================================================
// Perimeter Line Tracing (marching squares)
// ============================================================================
// The distance field is sampled at cell centres. For every square of four neighbouring centres, the level
// line crosses each side whose two corners lie on opposite sides of the level; the crossing point is found by
// linear interpolation. Each crossing is identified by the side it is on, so the pieces from neighbouring
// squares join up exactly and can be chained into closed loops.
//
// Orientation: walking round a square's corners in a fixed order, the sides where the line crosses alternate
// between "leaving the fire side" and "entering the fire side". Every piece is drawn from a leaving crossing
// to an entering crossing, which keeps the fire side on the same hand for every piece of every loop. The
// result is flipped once at the end if needed so the loops run anticlockwise on screen.
namespace {
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
}

std::vector<std::vector<Vec2>> tracePerimeterLine(const DistanceTransform& perimEdt, float level) {
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

    // Each sample is a corner of up to four squares. Rather than work f out (a square root) four times per
    // sample, keep the values of the square row's top and bottom sample rows and slide them down one row at a
    // time. The values are the same f(x, y) as before, so the traced line is exactly the same.
    // Samples run from x = win.x0 - 1 to win.x1 inclusive: entry k of a row buffer is sample x = win.x0 - 1 + k.
    const int sx0 = win.x0 - 1, nSamples = win.x1 - win.x0 + 2;
    std::vector<double> topRow(nSamples), bottomRow(nSamples);
    auto fillRow = [&](std::vector<double>& row, int y) {
        for (int k = 0; k < nSamples; ++k) row[k] = f(sx0 + k, y);
    };
    fillRow(bottomRow, win.y0 - 1);   // Becomes the top row of the first square row below

    for (int y = win.y0 - 1; y < win.y1; ++y) {
        // Slide down: the previous bottom row is this square row's top row
        std::swap(topRow, bottomRow);
        fillRow(bottomRow, y + 1);
        for (int x = win.x0 - 1; x < win.x1; ++x) {
            const int k = x - sx0;
            // Corner values in walking order (top-left, top-right, bottom-right, bottom-left), as SQ_DX / SQ_DY
            double v[4] = { topRow[k], topRow[k + 1], bottomRow[k + 1], bottomRow[k] };
            bool in[4];
            int nIn = 0;
            for (int c = 0; c < 4; ++c) {
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

// Shoelace formula over the closed polyline
double loopSignedArea(const std::vector<Vec2>& loop) {
    double a = 0.0;
    for (size_t i = 0, n = loop.size(); i < n; ++i) {
        const Vec2 p = loop[i], q = loop[(i + 1) % n];
        a += static_cast<double>(p.x) * q.y - static_cast<double>(q.x) * p.y;
    }
    return 0.5 * a;
}

float loopLength(const std::vector<Vec2>& loop) {
    float len = 0.0f;
    for (size_t i = 0, n = loop.size(); i < n; ++i) len += std::sqrt(sqDist(loop[i], loop[(i + 1) % n]));
    return len;
}
