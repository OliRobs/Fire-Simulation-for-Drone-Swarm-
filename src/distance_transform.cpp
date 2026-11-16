#include "distance_transform.h"
#include <algorithm>

// ============================================================================
// Euclidean Distance Transform (EDT)
//
// Computes the exact squared Euclidean Distance Transform of a 2D binary grid.
// Implements the linear-time O(N) separable algorithm by Felzenszwalb & Huttenlocher.
// ============================================================================

// Constructor: Initializes grid dimensions and pre-allocates workspace buffers
DistanceTransform::DistanceTransform(int w, int h)
    : w_(w), h_(h), sq_(static_cast<size_t>(w)* h) {
    const int m = std::max(w, h);
    // Allocate 1D buffers used by the 1D EDT pass (sized to max dimension m)
    f_.resize(m); d_.resize(m); v_.resize(m); z_.resize(m + 1);
}

// ----------------------------------------------------------------------------
// Compute 2D Distance Transform
// ----------------------------------------------------------------------------
// Inputs a binary mask (`sources`), where non-zero indicates a target/obstacle cell (distance = 0).
// Populates internal `sq_` vector with squared Euclidean distance to the nearest source cell.
void DistanceTransform::compute(const std::vector<std::uint8_t>& sources) {
    // --- Pass 1: distance to the nearest source in the same COLUMN ---
    // The input is binary (0 or INF), so the 1D transform down a column is just the squared gap to the
    // nearest source above or below. That is found with two row-by-row sweeps (down, then up), which read
    // memory in order: much faster than walking each column (a 2000-cell stride per step on a 2000 x 2000
    // grid), and the result is identical to running edt1d on every column.
    lastSource_.assign(w_, -1);   // Row of the last source seen in each column, -1 = none yet
    for (int y = 0; y < h_; ++y) {
        const int row = y * w_;
        for (int x = 0; x < w_; ++x) {
            if (sources[row + x]) lastSource_[x] = y;
            const int gap = y - lastSource_[x];
            sq_[row + x] = lastSource_[x] >= 0 ? double(gap) * gap : INF;
        }
    }
    lastSource_.assign(w_, -1);   // Now the nearest source BELOW, sweeping up from the bottom row
    for (int y = h_ - 1; y >= 0; --y) {
        const int row = y * w_;
        for (int x = 0; x < w_; ++x) {
            if (sources[row + x]) lastSource_[x] = y;
            if (lastSource_[x] >= 0) {
                const int gap = lastSource_[x] - y;
                sq_[row + x] = std::min(sq_[row + x], double(gap) * gap);
            }
        }
    }

    // --- Pass 2: Compute 1D Euclidean Distance Transform along rows ---
    // Combines the column distances: the lower envelope of parabolas gives the exact 2D distance.
    for (int y = 0; y < h_; ++y) {
        for (int x = 0; x < w_; ++x) {
            f_[x] = sq_[y * w_ + x]; // Load intermediate column EDT values
        }
        edt1d(w_); // Perform 1D EDT across row y
        for (int x = 0; x < w_; ++x) {
            sq_[y * w_ + x] = d_[x]; // Store final squared Euclidean distances
        }
    }
}

// ----------------------------------------------------------------------------
// 1D Euclidean Distance Transform (Felzenszwalb & Huttenlocher Algorithm)
// ----------------------------------------------------------------------------
// Calculates lower envelope of parabolas $f(p) + (q - p)^2$ in linear $O(n)$ time.
void DistanceTransform::edt1d(int n) {
    int k = 0;
    v_[0] = 0;        // Array of parabola locations forming the lower envelope
    z_[0] = -1e30;    // Array of intersection locations between adjacent parabolas
    z_[1] = 1e30;

    // Phase 1: Compute lower envelope boundary intersections
    for (int q = 1; q < n; ++q) {
        double s;
        for (;;) {
            const int p = v_[k];
            // Calculate intersection point 's' of parabolas starting at p and q
            s = ((f_[q] + double(q) * q) - (f_[p] + double(p) * p)) / (2.0 * q - 2.0 * p);

            // If intersection is left of previous boundary, top parabola is dominated; pop it
            if (s <= z_[k]) {
                --k;
            }
            else {
                break;
            }
        }
        // Push new parabola location and intersection boundary
        ++k;
        v_[k] = q;
        z_[k] = s;
        z_[k + 1] = 1e30;
    }

    // Phase 2: Fill distance values by evaluating the envelope at each grid index
    k = 0;
    for (int q = 0; q < n; ++q) {
        // Find which parabola in the lower envelope governs index q
        while (z_[k + 1] < q) {
            ++k;
        }
        const double diff = q - v_[k];
        // Calculate squared distance from q to nearest source vertex
        d_[q] = diff * diff + f_[v_[k]];
    }
}