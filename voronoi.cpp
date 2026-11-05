#include "voronoi.h"
#include <algorithm>
#include <limits>

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

    // ------------------------------------------------------------------------
    // Step 1: Farthest-Point Initial Seed Selection
    // ------------------------------------------------------------------------
    {
        // Tracks squared distance from each cell to its closest already-selected site
        std::vector<float> minD(cells.size(), std::numeric_limits<float>::max());

        // Pick the first seed centroid randomly from the set of valid cells
        sites.push_back(cellCentre(cells[rng() % cells.size()]));

        // Iteratively pick subsequent seeds at the location farthest from existing seeds
        while (static_cast<int>(sites.size()) < S) {
            size_t far = 0;
            float farD = -1.0f;

            for (size_t i = 0; i < cells.size(); ++i) {
                // Update minimum distance to any chosen site
                minD[i] = std::min(minD[i], sqDist(cellCentre(cells[i]), sites.back()));

                // Track cell with maximum distance from all current sites
                if (minD[i] > farD) {
                    farD = minD[i];
                    far = i;
                }
            }
            sites.push_back(cellCentre(cells[far]));
        }
    }

    // ------------------------------------------------------------------------
    // Step 2: Nearest-Site Voronoi Assignment Lambda
    // ------------------------------------------------------------------------
    // Assigns each cell to the index of its nearest site centroid (Euclidean distance)
    auto assignAll = [&]() {
        for (size_t i = 0; i < cells.size(); ++i) {
            const Vec2 p = cellCentre(cells[i]);
            int bestSite = 0;
            float bd = std::numeric_limits<float>::max();

            for (int s = 0; s < S; ++s) {
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
            sum[assign[i]] += cellCentre(cells[i]);
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