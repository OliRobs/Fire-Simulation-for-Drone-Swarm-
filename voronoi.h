#pragma once
#include <random>
#include <vector>
#include "types.h"

struct Voronoi {
    std::vector<Vec2> sites;
    std::vector<int> assign;   // Region index for each entry of the input cell list
};

// Splits the given cells (grid indices) into numSites Voronoi regions: farthest-point seeding followed by
// Lloyd relaxation so the regions come out roughly even.
Voronoi buildVoronoi(const std::vector<int>& cells, int numSites, std::mt19937& rng);
