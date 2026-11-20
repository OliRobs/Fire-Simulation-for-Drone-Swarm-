#include "distance_transform.h"
#include <algorithm>

// ============================================================================
// Euclidean Distance Transform (EDT)
//
// Computes the exact squared Euclidean Distance Transform of a 2D binary grid.
// Implements the linear-time O(N) separable algorithm by Felzenszwalb & Huttenlocher.
// ============================================================================

// Constructor: Initializes grid dimensions and pre-allocates workspace buffers. Every cell starts "very far"
// (INF), so with a range limit the cells never inside a window read correctly without being written.
DistanceTransform::DistanceTransform(int w, int h, int maxDistance)
    : w_(w), h_(h), maxDistance_(maxDistance), sq_(static_cast<size_t>(w) * h, INF) {
    const int m = std::max(w, h);
    // Allocate 1D buffers used by the 1D EDT pass (sized to max dimension m)
    f_.resize(m); d_.resize(m); v_.resize(m); z_.resize(m + 1);
}

// ----------------------------------------------------------------------------
// Compute 2D Distance Transform
// ----------------------------------------------------------------------------
// Inputs a binary mask (`sources`), where non-zero indicates a target/obstacle cell (distance = 0).
// Populates internal `sq_` vector with squared Euclidean distance to the nearest source cell.
// With a range limit (maxDistance_ > 0) only a window round the sources is worked out; see the header.
void DistanceTransform::compute(const std::vector<std::uint8_t>& sources) {
    if (maxDistance_ <= 0) {   // No range limit: the bounding box is not needed, the whole grid is worked out
        computeIn(sources, GridRect{ 0, 0, w_, h_ });
        return;
    }
    // --- Pass 0: bounding box of the sources ---
    // A row with no source is skipped with one std::find_if over its bytes. Callers that already know the
    // box (the fire simulation tracks it) use the other overload and skip this 4 MB scan.
    GridRect box;
    for (int y = 0; y < h_; ++y) {
        const std::uint8_t* row = sources.data() + static_cast<size_t>(y) * w_;
        const std::uint8_t* first = std::find_if(row, row + w_, [](std::uint8_t v) { return v != 0; });
        if (first == row + w_) continue;
        const std::uint8_t* last = row + w_ - 1;
        while (*last == 0) --last;   // Stops at 'first' at the latest
        box.include(static_cast<int>(first - row), y);
        box.include(static_cast<int>(last - row), y);
    }
    computeIn(sources, box);
}

// ----------------------------------------------------------------------------
// Compute with a known source rectangle
// ----------------------------------------------------------------------------
// The caller promises every source lies inside sourceBounds (see the header).
void DistanceTransform::compute(const std::vector<std::uint8_t>& sources, const GridRect& sourceBounds) {
    // Clip to the grid, in case the caller's rectangle runs off it
    const GridRect clipped{ std::max(0, sourceBounds.x0), std::max(0, sourceBounds.y0),
                            std::min(w_, sourceBounds.x1), std::min(h_, sourceBounds.y1) };
    computeIn(sources, maxDistance_ > 0 ? clipped : GridRect{ 0, 0, w_, h_ });
}

// ----------------------------------------------------------------------------
// Shared body of both compute() overloads
// ----------------------------------------------------------------------------
// sourceBounds: a rectangle holding every source cell (the whole grid when there is no range limit).
void DistanceTransform::computeIn(const std::vector<std::uint8_t>& sources, const GridRect& sourceBounds) {
    ++version_;   // Every compute() goes through here: readers can tell the field has been recomputed
    // --- Choose the window of cells to work out ---
    GridRect win{ 0, 0, w_, h_ };   // No range limit: the whole grid
    if (maxDistance_ > 0) {
        // Grow the box by the range (+1 so a threshold of exactly maxDistance is still exact) and clip to the grid
        const int r = maxDistance_ + 1;
        win = sourceBounds.empty() ? GridRect{}   // No sources at all: nothing to work out, everything is INF
                                   : GridRect{ std::max(0, sourceBounds.x0 - r), std::max(0, sourceBounds.y0 - r),
                                               std::min(w_, sourceBounds.x1 + r), std::min(h_, sourceBounds.y1 + r) };

        // Cells the previous window wrote but this one may not: back to INF. (The fire only grows, so the new
        // window normally contains the old one and this just rewrites cells that are about to be overwritten.)
        for (int y = last_.y0; y < last_.y1; ++y) {
            std::fill(sq_.begin() + static_cast<size_t>(y) * w_ + last_.x0,
                      sq_.begin() + static_cast<size_t>(y) * w_ + last_.x1, INF);
        }
        last_ = win;
    }
    const int wx0 = win.x0, wx1 = win.x1, wy0 = win.y0, wy1 = win.y1;
    if (wx1 <= wx0 || wy1 <= wy0) return;
    // All sources lie inside the window, so working only on its rows and columns gives the same exact distances
    // as the whole-grid transform for every cell in it.

    // --- Pass 1: distance to the nearest source in the same COLUMN ---
    // The input is binary (0 or INF), so the 1D transform down a column is just the squared gap to the
    // nearest source above or below. That is found with two row-by-row sweeps (down, then up), which read
    // memory in order: much faster than walking each column (a 2000-cell stride per step on a 2000 x 2000
    // grid), and the result is identical to running edt1d on every column.
    lastSource_.assign(w_, -1);   // Row of the last source seen in each column, -1 = none yet
    for (int y = wy0; y < wy1; ++y) {
        const int row = y * w_;
        for (int x = wx0; x < wx1; ++x) {
            if (sources[row + x]) lastSource_[x] = y;
            const int gap = y - lastSource_[x];
            sq_[row + x] = lastSource_[x] >= 0 ? double(gap) * gap : INF;
        }
    }
    lastSource_.assign(w_, -1);   // Now the nearest source BELOW, sweeping up from the bottom row
    for (int y = wy1 - 1; y >= wy0; --y) {
        const int row = y * w_;
        for (int x = wx0; x < wx1; ++x) {
            if (sources[row + x]) lastSource_[x] = y;
            if (lastSource_[x] >= 0) {
                const int gap = lastSource_[x] - y;
                sq_[row + x] = std::min(sq_[row + x], double(gap) * gap);
            }
        }
    }

    // --- Pass 2: Compute 1D Euclidean Distance Transform along rows ---
    // Combines the column distances: the lower envelope of parabolas gives the exact 2D distance.
    // Only the window's columns: outside them the column distances are all INF, which never wins.
    const int n = wx1 - wx0;
    for (int y = wy0; y < wy1; ++y) {
        const int row = y * w_ + wx0;
        for (int x = 0; x < n; ++x) {
            f_[x] = sq_[row + x]; // Load intermediate column EDT values
        }
        edt1d(n); // Perform 1D EDT across the window's part of row y
        for (int x = 0; x < n; ++x) {
            sq_[row + x] = d_[x]; // Store final squared Euclidean distances
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