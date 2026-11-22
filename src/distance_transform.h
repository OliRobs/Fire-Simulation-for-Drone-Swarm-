#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>
#include "types.h"

// Exact Euclidean distance transform (Felzenszwalb & Huttenlocher), O(N).
// Gives each cell its squared distance to the nearest source cell.
//
// Optional range limit: when maxDistance > 0, only the cells within maxDistance of the sources' bounding box are
// worked out (the window). Distances inside the window are still exact; cells outside it, which are all more
// than maxDistance from every source, read as "very far" (INF). The fire covers a tiny part of the 2000 x 2000
// grid for most of a run, so this cuts the per-frame cost from the whole grid to a window round the fire.
// Use it only where every reader compares the distance against thresholds no bigger than maxDistance.
//
// Incremental mode (for the live fire, whose cells only ever switch on): distances are exact out to
// maxDistance + 1 cells and every cell further away reads INF ("very far"), instead of its exact distance. Every
// reader of the live field compares against thresholds within maxDistance, so it sees exactly the same answers.
// In exchange, compute() only has to deal with the source cells that switched on since the previous call: each
// lowers the distances inside a disc of radius maxDistance + 1 round it, and nothing else can change. That is a few
// thousand cells per new fire cell instead of three passes over the whole window. It falls back to a full
// transform (then clamped) on the first call, when a source cell has switched OFF, and when so many cells switched
// on at once that the full transform is cheaper.
class DistanceTransform {
public:
    enum class Mode { Exact, IncrementalClamped };

    // maxDistance: cells; 0 = no limit (the whole grid, the original behaviour). IncrementalClamped needs a limit.
    DistanceTransform(int w, int h, int maxDistance = 0, Mode mode = Mode::Exact);

    // sources[i] != 0 marks a source cell. Read results via squaredDistance().
    void compute(const std::vector<std::uint8_t>& sources);

    // Same, when the caller already knows a rectangle that contains every source (e.g. FireSimulation's
    // fireBounds(), which also contains the perimeter's committed fire). With a range limit this skips the
    // whole-grid scan for the sources' bounding box, the main per-frame cost once the window is small. The
    // result is exact as long as no source lies outside sourceBounds; a looser rectangle only makes the
    // window bigger. Without a range limit sourceBounds is not needed and is ignored.
    void compute(const std::vector<std::uint8_t>& sources, const GridRect& sourceBounds);

    double squaredDistance(int idx) const { return sq_[idx]; }

    // A lower bound on the (not squared) distance from cell idx to the nearest source cell: the exact distance
    // where squaredDistance() holds it, and where it reads INF ("very far") the least distance that can mean. With
    // a range limit INF means more than maxDistance + 1 cells from every source (outside the window, or beyond the
    // incremental clamp); without one it only appears when there are no sources at all.
    // Because the true distance changes by at most d between two cells d apart (triangle inequality), a reader can
    // use this to know that every cell within (bound - threshold) of idx is also beyond a threshold without
    // reading them (SafePlanner::lineClear does).
    double distanceLowerBound(int idx) const {
        const double far = maxDistance_ > 0 ? maxDistance_ + 1.0 : INF;
        return std::min(std::sqrt(sq_[idx]), far);
    }

    // Cells the last compute() worked out. Every cell outside it reads INF, i.e. is more than maxDistance from
    // every source, so a reader looking for cells nearer than that can loop over this rectangle only. The
    // whole grid when there is no range limit.
    GridRect window() const { return maxDistance_ > 0 ? last_ : GridRect{ 0, 0, w_, h_ }; }

    // Number of times compute() has run. A reader that keeps something worked out from these distances can
    // compare versions to tell whether the field may have changed since (e.g. the mission controller keeps the
    // ACO tour until the Perimeter, which only recomputes its transform when a section moves, changes).
    std::uint64_t version() const { return version_; }

private:
    static constexpr double INF = 1e9; // >> max possible squared distance on the grid
    void edt1d(int n);

    // Shared by both compute() overloads: works out the window grown from the sources' bounding box
    void computeIn(const std::vector<std::uint8_t>& sources, const GridRect& sourceBounds);
    // Incremental mode: update from the source cells that switched on, or fall back to computeIn + clamp
    void computeIncremental(const std::vector<std::uint8_t>& sources, const GridRect& sourceBounds);
    // Incremental mode: lowers every distance within the clamp radius of source cell (sx, sy)
    void addSource(int sx, int sy);

    int w_, h_;
    int maxDistance_;               // 0 = whole grid
    std::uint64_t version_ = 0;     // compute() calls so far (see version())
    GridRect last_;                 // Window [x0, x1) x [y0, y1) written by the previous compute(), reset to INF before the next one
    std::vector<double> sq_, f_, d_, z_;
    std::vector<int> v_;
    std::vector<int> lastSource_;   // Pass 1 scratch: row of the nearest source so far in each column

    // --- Incremental mode only ---
    bool incremental_ = false;
    bool haveSeen_ = false;              // seen_ holds the sources of the previous compute()
    std::vector<std::uint8_t> seen_;     // Copy of the sources at the previous compute() (to find new ones)
    GridRect seenBox_;                   // Every non-zero cell of seen_ lies inside this rectangle
    std::vector<int> newSources_;        // Scratch: cells that switched on since the previous compute()
};
