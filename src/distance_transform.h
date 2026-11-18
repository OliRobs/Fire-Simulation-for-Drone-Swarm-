#pragma once
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
class DistanceTransform {
public:
    // maxDistance: cells; 0 = no limit (the whole grid, the original behaviour)
    DistanceTransform(int w, int h, int maxDistance = 0);

    // sources[i] != 0 marks a source cell. Read results via squaredDistance().
    void compute(const std::vector<std::uint8_t>& sources);

    // Same, when the caller already knows a rectangle that contains every source (e.g. FireSimulation's
    // fireBounds(), which also contains the perimeter's committed fire). With a range limit this skips the
    // whole-grid scan for the sources' bounding box, the main per-frame cost once the window is small. The
    // result is exact as long as no source lies outside sourceBounds; a looser rectangle only makes the
    // window bigger. Without a range limit sourceBounds is not needed and is ignored.
    void compute(const std::vector<std::uint8_t>& sources, const GridRect& sourceBounds);

    double squaredDistance(int idx) const { return sq_[idx]; }

    // Cells the last compute() worked out. Every cell outside it reads INF, i.e. is more than maxDistance from
    // every source, so a reader looking for cells nearer than that can loop over this rectangle only. The
    // whole grid when there is no range limit.
    GridRect window() const { return maxDistance_ > 0 ? last_ : GridRect{ 0, 0, w_, h_ }; }

private:
    static constexpr double INF = 1e9; // >> max possible squared distance on the grid
    void edt1d(int n);

    // Shared by both compute() overloads: works out the window grown from the sources' bounding box
    void computeIn(const std::vector<std::uint8_t>& sources, const GridRect& sourceBounds);

    int w_, h_;
    int maxDistance_;               // 0 = whole grid
    GridRect last_;                 // Window [x0, x1) x [y0, y1) written by the previous compute(), reset to INF before the next one
    std::vector<double> sq_, f_, d_, z_;
    std::vector<int> v_;
    std::vector<int> lastSource_;   // Pass 1 scratch: row of the nearest source so far in each column
};
