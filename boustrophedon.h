#pragma once
#include <cstdint>
#include <vector>
#include "mission_types.h"

// Boustrophedon ("as the ox ploughs") sweep of one Voronoi region's inner-zone cells: parallel lanes spaced
// one sensor swath apart, direction alternating each lane. Appends the waypoints to route.
void appendBoustrophedon(int r, const Region& rg, const std::vector<int>& label,
                         const std::vector<std::uint8_t>& zone, std::vector<Vec2>& route);
