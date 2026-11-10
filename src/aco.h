#pragma once
#include <random>
#include <vector>
#include "types.h"

// Ant colony optimisation for the closed tour (TSP) through a set of points.
// Returns the visiting order as indices into pts.
std::vector<int> solveACO(const std::vector<Vec2>& pts, std::mt19937& rng);
