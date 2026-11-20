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

    // Cell centres, worked out once. The steps below visit cells many times (seeding, then once per Lloyd
    // pass), and cellCentre() does an integer divide and modulo per call.
    std::vector<Vec2> pts(cells.size());
    for (size_t i = 0; i < cells.size(); ++i) pts[i] = cellCentre(cells[i]);

    // ------------------------------------------------------------------------
    // Blocks of cells, shared by the seeding (Step 1) and the assignment (Step 2)
    // ------------------------------------------------------------------------
    // The cells' bounding box is cut into BUCKET x BUCKET blocks. Both steps use a block's rectangle to rule
    // out, with one cheap bound per block, work that cannot change any cell inside it.
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

    // Cells of each block, in increasing cell-list order: block b holds blockCells[blockStart[b] ..
    // blockStart[b + 1]). Built with a counting sort, so the order inside each block is the input order.
    std::vector<int> blockStart(static_cast<size_t>(bw) * bh + 1, 0);
    for (size_t i = 0; i < cells.size(); ++i) ++blockStart[bucketOf[i] + 1];
    for (int b = 0; b < bw * bh; ++b) blockStart[b + 1] += blockStart[b];
    std::vector<int> blockCells(cells.size());
    {
        std::vector<int> fill(blockStart.begin(), blockStart.end() - 1);
        for (size_t i = 0; i < cells.size(); ++i) blockCells[fill[bucketOf[i]]++] = static_cast<int>(i);
    }

    // Rectangle spanned by the cell centres of block b (centres are at +0.5). Every cell centre in the block
    // lies inside [rx0, rx1] x [ry0, ry1].
    auto blockRect = [&](int b, float& rx0, float& rx1, float& ry0, float& ry1) {
        rx0 = minX + (b % bw) * BUCKET + 0.5f; rx1 = rx0 + (BUCKET - 1);
        ry0 = minY + (b / bw) * BUCKET + 0.5f; ry1 = ry0 + (BUCKET - 1);
        };

    // Smallest squared distance from point q to the rectangle (0 when q is inside it): a lower bound on
    // sqDist(q, p) for every cell centre p in the block
    auto rectLowerBound = [](const Vec2& q, float rx0, float rx1, float ry0, float ry1) {
        const float nx = q.x < rx0 ? rx0 - q.x : (q.x > rx1 ? q.x - rx1 : 0.0f);
        const float ny = q.y < ry0 ? ry0 - q.y : (q.y > ry1 ? q.y - ry1 : 0.0f);
        return nx * nx + ny * ny;
        };

    // ------------------------------------------------------------------------
    // Step 1: Farthest-Point Initial Seed Selection (bucketed, exact)
    // ------------------------------------------------------------------------
    // Each new seed is the cell farthest from all seeds chosen so far (the first such cell in list order on a
    // tie). The plain method updates every cell's nearest-seed distance once per seed: cells x S distance
    // tests, the largest part of a re-plan once the assignment was bucketed.
    //
    // A new seed q only lowers the distance of cells that are nearer to q than to every earlier seed. Each
    // block keeps the largest nearest-seed distance of its cells (blockMax); if q is further from the whole
    // block rectangle than that, no cell in the block can change, so the block is skipped. Blocks that are
    // updated also re-find their farthest cell. The farthest cell overall is then the block maximum with the
    // smallest cell-list index on a tie, which is exactly the cell the plain scan picks: the same float
    // distances are compared in the same way, only cells that cannot change are left out.
    {
        // Squared distance from each cell to its closest already-selected site
        std::vector<float> minD(cells.size(), std::numeric_limits<float>::max());
        // Per block: largest minD of its cells, and the first cell (cell-list index) that has it
        std::vector<float> blockMax(static_cast<size_t>(bw) * bh, -1.0f);
        std::vector<int> blockArg(static_cast<size_t>(bw) * bh, -1);
        for (int b = 0; b < bw * bh; ++b) {
            if (blockStart[b] == blockStart[b + 1]) continue;   // Empty block: stays at -1, never picked
            blockMax[b] = std::numeric_limits<float>::max();
            blockArg[b] = blockCells[blockStart[b]];
        }

        // Pick the first seed centroid randomly from the set of valid cells
        sites.push_back(pts[rng() % cells.size()]);

        // Iteratively pick subsequent seeds at the location farthest from existing seeds
        while (static_cast<int>(sites.size()) < S) {
            const Vec2 q = sites.back();
            int far = -1;          // Cell-list index of the farthest cell so far
            float farD = -1.0f;

            for (int b = 0; b < bw * bh; ++b) {
                if (blockArg[b] < 0) continue;   // Empty block
                float rx0, rx1, ry0, ry1;
                blockRect(b, rx0, rx1, ry0, ry1);
                // Skip the block when every cell in it is certainly nearer an earlier seed than q. The margin
                // (lower bound shrunk a little) absorbs float rounding, so a block is only skipped when the
                // exact distance test would leave every one of its cells unchanged.
                const float lb = rectLowerBound(q, rx0, rx1, ry0, ry1);
                if (lb * (1.0f - 1e-4f) - 1e-3f <= blockMax[b]) {
                    // Update the block's cells, in list order, and re-find its farthest cell
                    float bm = -1.0f;
                    int ba = -1;
                    for (int c = blockStart[b], e = blockStart[b + 1]; c < e; ++c) {
                        const int i = blockCells[c];
                        minD[i] = std::min(minD[i], sqDist(pts[i], q));
                        if (minD[i] > bm) { bm = minD[i]; ba = i; }
                    }
                    blockMax[b] = bm;
                    blockArg[b] = ba;
                }
                // Overall farthest cell: largest distance, then the smallest cell-list index
                if (blockMax[b] > farD || (blockMax[b] == farD && blockArg[b] < far)) {
                    farD = blockMax[b];
                    far = blockArg[b];
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
    std::vector<int> candStart(static_cast<size_t>(bw) * bh + 1);   // Block b's candidates: [candStart[b], candStart[b + 1])
    std::vector<int> cand;                                          // Candidate site lists, all blocks back to back
    std::vector<float> lower(S), upper(S);

    auto assignAll = [&]() {
        // --- 2a: candidate sites per block, from the block's rectangle of cell centres ---
        cand.clear();
        for (int b = 0; b < bw * bh; ++b) {
            candStart[b] = static_cast<int>(cand.size());
            if (!used[b]) continue;
            // Rectangle spanned by the cell centres of block b
            float rx0, rx1, ry0, ry1;
            blockRect(b, rx0, rx1, ry0, ry1);
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
        // Block by block, so a block's candidate list stays in cache while its cells are done. Each cell's answer
        // depends only on its own position and the sites, not on the order cells are visited in, so this gives
        // exactly the assignment of a pass in cell-list order.
        for (int b = 0; b < bw * bh; ++b) {
            const int c0 = candStart[b], c1 = candStart[b + 1];
            if (c0 == c1) continue;   // Empty block (no cells, no candidates)
            // Only one site can be nearest anywhere in this block: every cell in it goes to that site, with no
            // distance tests (a full scan would pick it too, as the only site not beaten everywhere in the block)
            if (c1 - c0 == 1) {
                for (int k = blockStart[b], e = blockStart[b + 1]; k < e; ++k) assign[blockCells[k]] = cand[c0];
                continue;
            }
            for (int k = blockStart[b], e = blockStart[b + 1]; k < e; ++k) {
                const int i = blockCells[k];
                const Vec2 p = pts[i];
                int bestSite = 0;
                float bd = std::numeric_limits<float>::max();
                for (int c = c0; c < c1; ++c) {
                    const int s = cand[c];
                    const float dd = sqDist(p, sites[s]);
                    if (dd < bd) {
                        bd = dd;
                        bestSite = s;
                    }
                }
                assign[i] = bestSite;
            }
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