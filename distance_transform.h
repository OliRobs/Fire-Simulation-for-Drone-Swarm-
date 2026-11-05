#pragma once
#include <cstdint>
#include <vector>

// Exact Euclidean distance transform (Felzenszwalb & Huttenlocher), O(N).
// Gives each cell its squared distance to the nearest source cell.
class DistanceTransform {
public:
    DistanceTransform(int w, int h);

    // sources[i] != 0 marks a source cell. Read results via squaredDistance().
    void compute(const std::vector<std::uint8_t>& sources);
    double squaredDistance(int idx) const { return sq_[idx]; }

private:
    static constexpr double INF = 1e9; // >> max possible squared distance on the grid
    void edt1d(int n);

    int w_, h_;
    std::vector<double> sq_, f_, d_, z_;
    std::vector<int> v_;
};
