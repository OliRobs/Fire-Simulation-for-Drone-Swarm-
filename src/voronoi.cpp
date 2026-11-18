#include "voronoi.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstddef>

// ============================================================================
// Voronoi Spatial Partitioning & Lloyd's Relaxation Engine
//
// Splits a set of grid cells into `numSites` spatially balanced regions:
// 1. Farthest-Point Seeding (k-means++ style initialization): Selects initial 
//    seed centroids spaced as far apart as possible to prevent cluster overlap.
// 2. Voronoi Cell Assignment: Maps each grid cell to its nearest centroid.
// 3. Lloyd's Relaxation: Iteratively moves site centroids to the center of 
//    mass of their assigned cell regions to create uniform, well-shaped sectors.
// ============================================================================
Voronoi buildVoronoi(const std::vector<int>& cells, int numSites, std::mt19937& rng) {
    const int S = numSites;
    Voronoi v;
    std::vector<Vec2>& sites = v.sites;
    std::vector<int>& assign = v.assign;
    assign.assign(cells.size(), 0);
    if (cells.empty() || S <= 0) return v;

    // Cell centres, worked out once. Every step below loops over all cells many times (once per seed, then
    // once per Lloyd pass), and cellCentre() does an integer divide and modulo per call.
    std::vector<Vec2> pts(cells.size());
    for (size_t i = 0; i < cells.size(); ++i) pts[i] = cellCentre(cells[i]);

    // ------------------------------------------------------------------------
    // Step 1: Farthest-Point Initial Seed Selection
    // ------------------------------------------------------------------------
    {
        // Tracks squared distance from each cell to its closest already-selected site
        std::vector<float> minD(cells.size(), std::numeric_limits<float>::max());

        // Pick the first seed centroid randomly from the set of valid cells
        sites.push_back(pts[rng() % cells.size()]);

        // Iteratively pick subsequent seeds at the location farthest from existing seeds
        while (static_cast<int>(sites.size()) < S) {
            size_t far = 0;
            float farD = -1.0f;

            for (size_t i = 0; i < cells.size(); ++i) {
                // Update minimum distance to any chosen site
                minD[i] = std::min(minD[i], sqDist(pts[i], sites.back()));

                // Track cell with maximum distance from all current sites
                if (minD[i] > farD) {
                    farD = minD[i];
                    far = i;
                }
            }
            sites.push_back(pts[far]);
        }
    }

    // ------------------------------------------------------------------------
    // Step 2: Nearest-Site Voronoi Assignment (bucketed, exact)
    // ------------------------------------------------------------------------
    // Assigns each cell to the index of its nearest site (Euclidean distance; on a tie, the lowest site
    // index). Checking all S sites for every cell costs cells x S distance tests per pass (about 20 M for a
    // large fire with 64 sites, five passes per plan), which was most of a re-plan.
    //
    // Instead the cells' bounding box is cut into BUCKET x BUCKET blocks. For each block we find, from the
    // block's rectangle alone, the sites that could be nearest to some point in it: the nearest any site can
    // be is its distance to the rectangle (lower bound), the furthest is its distance to the far corner
    // (upper bound). A site whose lower bound exceeds the smallest upper bound of all sites is beaten
    // everywhere in the block by that other site, so it can never win. Only the remaining candidates (a
    // handful near the block) are tested per cell, in increasing site order with the same comparison as a
    // full scan, so the answer, ties included, is exactly the same as testing every site.
    constexpr int BUCKET = 16;   // Block side in cells: small enough to prune well, big enough to be cheap
    int minX = GRID_WIDTH, minY = GRID_HEIGHT, maxX = -1, maxY = -1;
    for (int idx : cells) {
        const int x = idx % GRID_WIDTH, y = idx / GRID_WIDTH;
        minX = std::min(minX, x); maxX = std::max(maxX, x);
        minY = std::min(minY, y); maxY = std::max(maxY, y);
    }
    const int bw = (maxX - minX) / BUCKET + 1, bh = (maxY - minY) / BUCKET + 1;

    // Block of each cell (fixed: the cells don't move, only the sites do)
    std::vector<int> bucketOf(cells.size());
    std::vector<char> used(static_cast<size_t>(bw) * bh, 0);   // Blocks that hold at least one cell
    for (size_t i = 0; i < cells.size(); ++i) {
        const int x = cells[i] % GRID_WIDTH, y = cells[i] / GRID_WIDTH;
        bucketOf[i] = ((y - minY) / BUCKET) * bw + (x - minX) / BUCKET;
        used[bucketOf[i]] = 1;
    }

    std::vector<int> candStart(static_cast<size_t>(bw) * bh + 1);   // Block b's candidates: [candStart[b], candStart[b + 1])
    std::vector<int> cand;                                          // Candidate site lists, all blocks back to back
    std::vector<float> lower(S), upper(S);

    auto assignAll = [&]() {
        // --- 2a: candidate sites per block, from the block's rectangle of cell centres ---
        cand.clear();
        for (int b = 0; b < bw * bh; ++b) {
            candStart[b] = static_cast<int>(cand.size());
            if (!used[b]) continue;
            // Rectangle spanned by the cell centres of block b (centres are at +0.5)
            const float rx0 = minX + (b % bw) * BUCKET + 0.5f, rx1 = rx0 + (BUCKET - 1);
            const float ry0 = minY + (b / bw) * BUCKET + 0.5f, ry1 = ry0 + (BUCKET - 1);
            float bestUpper = std::numeric_limits<float>::max();
            for (int s = 0; s < S; ++s) {
                const Vec2 q = sites[s];
                // Nearest point of the rectangle to the site (0 along an axis where the site is inside it)
                const float nx = q.x < rx0 ? rx0 - q.x : (q.x > rx1 ? q.x - rx1 : 0.0f);
                const float ny = q.y < ry0 ? ry0 - q.y : (q.y > ry1 ? q.y - ry1 : 0.0f);
                // Furthest corner of the rectangle from the site
                const float fx = std::max(std::fabs(q.x - rx0), std::fabs(q.x - rx1));
                const float fy = std::max(std::fabs(q.y - ry0), std::fabs(q.y - ry1));
                lower[s] = nx * nx + ny * ny;
                upper[s] = fx * fx + fy * fy;
                bestUpper = std::min(bestUpper, upper[s]);
            }
            // Keep every site that might win. The small relative margin absorbs float rounding, so a site
            // is only dropped when it is clearly further than another site for every cell in the block.
            const float keep = bestUpper * (1.0f + 1e-4f) + 1e-3f;
            for (int s = 0; s < S; ++s) {
                if (lower[s] <= keep) cand.push_back(s);
            }
        }
        candStart[bw * bh] = static_cast<int>(cand.size());

        // --- 2b: each cell against its block's candidates only (in increasing site order) ---
        for (size_t i = 0; i < cells.size(); ++i) {
            const Vec2 p = pts[i];
            int bestSite = 0;
            float bd = std::numeric_limits<float>::max();
            for (int c = candStart[bucketOf[i]], e = candStart[bucketOf[i] + 1]; c < e; ++c) {
                const int s = cand[c];
                const float dd = sqDist(p, sites[s]);
                if (dd < bd) {
                    bd = dd;
                    bestSite = s;
                }
            }
            assign[i] = bestSite;
        }
        };

    // ------------------------------------------------------------------------
    // Step 3: Lloyd's Relaxation Iterations (Centroid Refinement)
    // ------------------------------------------------------------------------
    // Iteratively reassigns cells and recomputes region centroids to smooth boundaries
    for (int it = 0; it < LLOYD_ITERATIONS; ++it) {
        assignAll(); // Assign cells to current site positions

        std::vector<Vec2> sum(S, Vec2(0.0f, 0.0f));
        std::vector<int> cnt(S, 0);

        // Accumulate center positions of all cells belonging to each site
        for (size_t i = 0; i < cells.size(); ++i) {
            sum[assign[i]] += pts[i];
            ++cnt[assign[i]];
        }

        // Relocate each site to the center-of-mass (mean position) of its assigned cells
        for (int s = 0; s < S; ++s) {
            if (cnt[s] > 0) {
                sites[s] = sum[s] / static_cast<float>(cnt[s]);
            }
        }
    }

    // Final assignment pass to lock in the refined Voronoi partition boundaries
    assignAll();
    return v;
}