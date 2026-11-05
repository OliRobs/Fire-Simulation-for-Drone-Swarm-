#include "aco.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

// ============================================================================
// Ant Colony Optimization (ACO) for the Traveling Salesperson Problem (TSP)
// 
// Algorithm Variant: Elitist Ant System (EAS)
// - Construction: Ants construct tours probabilistically based on edge pheromones
//   and heuristic info (1 / distance^3).
// - Evaporation: Pheromones decay uniformly across all edges per iteration.
// - Deposit: Every ant deposits pheromone (1 / tour_length) along its paths.
// - Elite Deposit: The global best tour receives additional pheromone reinforcement.
// ============================================================================
std::vector<int> solveACO(const std::vector<Vec2>& pts, std::mt19937& rng) {
    const int n = static_cast<int>(pts.size());
    std::vector<int> best(n);

    // Initialize default tour as identity permutation: [0, 1, 2, ..., n-1]
    std::iota(best.begin(), best.end(), 0);

    // Trivial cases: Tours with 3 or fewer points have fixed/equivalent path lengths
    if (n <= 3) return best;

    // ------------------------------------------------------------------------
    // Step 1: Precompute Distance Matrix & Heuristic Values (1 / dist^3)
    // ------------------------------------------------------------------------
    std::vector<float> d(static_cast<size_t>(n) * n), heur(static_cast<size_t>(n) * n, 0.0f);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            d[i * n + j] = std::sqrt(sqDist(pts[i], pts[j]));
            if (i != j) {
                // Eta (heuristic desirability): Inverse distance cubed heavily favors closer points
                const float h = 1.0f / (d[i * n + j] + 1e-3f);
                heur[i * n + j] = h * h * h;
            }
        }
    }

    // Helper lambda to compute total round-trip tour distance
    auto tourLength = [&](const std::vector<int>& t) {
        float len = 0.0f;
        for (int i = 0; i < n; ++i) len += d[t[i] * n + t[(i + 1) % n]];
        return len;
        };

    // ------------------------------------------------------------------------
    // Step 2: Seed Global Best with Nearest-Neighbor Greedy Heuristic
    // ------------------------------------------------------------------------
    {
        std::vector<char> seen(n, 0);
        int cur = 0;
        best[0] = 0; seen[0] = 1;
        for (int step = 1; step < n; ++step) {
            int next = -1; float nd = std::numeric_limits<float>::max();
            for (int j = 0; j < n; ++j)
                if (!seen[j] && d[cur * n + j] < nd) { nd = d[cur * n + j]; next = j; }
            best[step] = next; seen[next] = 1; cur = next;
        }
    }
    float bestLen = tourLength(best);

    // ------------------------------------------------------------------------
    // Step 3: Initialize ACO Data Structures & Initial Pheromone Trail (Tau)
    // ------------------------------------------------------------------------
    // Initial pheromone scale derived from greedy tour length to keep weights balanced
    std::vector<float> tau(static_cast<size_t>(n) * n, ACO_ANTS / std::max(bestLen, 1e-3f));
    std::vector<std::vector<int>> tours(ACO_ANTS, std::vector<int>(n));
    std::vector<float> lens(ACO_ANTS), w(n);
    std::vector<char> visited(n);
    std::uniform_real_distribution<float> u01(0.0f, 1.0f);

    // ------------------------------------------------------------------------
    // Step 4: Main Iteration Loop
    // ------------------------------------------------------------------------
    for (int it = 0; it < ACO_ITERATIONS; ++it) {

        // --- 4a. Construct Tours for All Ants ---
        for (int a = 0; a < ACO_ANTS; ++a) {
            std::fill(visited.begin(), visited.end(), 0);

            // Random starting location for ant 'a'
            int cur = static_cast<int>(rng() % n);
            tours[a][0] = cur; visited[cur] = 1;

            // Step through and choose remaining cities
            for (int step = 1; step < n; ++step) {
                float sum = 0.0f;

                // Calculate transition probability weights: tau^alpha * eta^beta (here alpha=1, beta=1 implicitly embedded in heur)
                for (int j = 0; j < n; ++j) {
                    w[j] = visited[j] ? 0.0f : tau[cur * n + j] * heur[cur * n + j];
                    sum += w[j];
                }

                int next = -1;
                // Roulette-Wheel Selection (Stochastic Proportional Selection)
                if (sum > 0.0f) {
                    float pick = u01(rng) * sum;
                    for (int j = 0; j < n; ++j) {
                        if (w[j] <= 0.0f) continue;
                        next = j; pick -= w[j];
                        if (pick <= 0.0f) break;
                    }
                }

                // Fallback deterministic selection in case of numerical underflow
                if (next < 0)
                    for (int j = 0; j < n; ++j) if (!visited[j]) { next = j; break; }

                tours[a][step] = next; visited[next] = 1; cur = next;
            }

            // Evaluate ant's constructed tour and update global best if improved
            lens[a] = tourLength(tours[a]);
            if (lens[a] < bestLen) { bestLen = lens[a]; best = tours[a]; }
        }

        // --- 4b. Pheromone Evaporation ---
        for (float& t : tau) t *= (1.0f - ACO_EVAPORATION);

        // Helper lambda to deposit pheromones symmetrically across tour edges
        auto deposit = [&](const std::vector<int>& t, float amount) {
            for (int i = 0; i < n; ++i) {
                const int a = t[i], b = t[(i + 1) % n];
                tau[a * n + b] += amount;
                tau[b * n + a] += amount;
            }
            };

        // --- 4c. Pheromone Updates (Standard Ants + Elite Reinforcement) ---
        for (int a = 0; a < ACO_ANTS; ++a) deposit(tours[a], 1.0f / lens[a]);
        deposit(best, ACO_ELITE_WEIGHT / bestLen);
    }

    return best; // Return best tour permutation found
}