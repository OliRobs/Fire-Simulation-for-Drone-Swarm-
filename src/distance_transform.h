#pragma once
#include <cstdint>
#include <vector>

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
    double squaredDistance(int idx) const { return sq_[idx]; }

private:
    static constexpr double INF = 1e9; // >> max possible squared distance on the grid
    void edt1d(int n);

    // Rectangle of cells [x0, x1) x [y0, y1) that the last compute() wrote
    struct Window { int x0 = 0, y0 = 0, x1 = 0, y1 = 0; };

    int w_, h_;
    int maxDistance_;               // 0 = whole grid
    Window last_;                   // Window written by the previous compute(), reset to INF before the next one
    std::vector<double> sq_, f_, d_, z_;
    std::vector<int> v_;
    std::vector<int> lastSource_;   // Pass 1 scratch: row of the nearest source so far in each column
};
