#pragma once
#include <limits>
#include <vector>
#include "types.h"

struct Region {
    int innerCount = 0, outerCount = 0;   // Cells inside / outside the perimeter boundary
    int minX = std::numeric_limits<int>::max(), maxX = -1; // Bounding box of the inner part
    int minY = std::numeric_limits<int>::max(), maxY = -1;
    Vec2 innerSum{ 0.0f, 0.0f }, outerSum{ 0.0f, 0.0f };
    Vec2 innerCentroid{ 0.0f, 0.0f }, outerCentroid{ 0.0f, 0.0f };
    Vec2 poi{ 0.0f, 0.0f };               // Outer-zone point of interest (visited by the ACO tour)
    float poiDist = std::numeric_limits<float>::max();
};

struct MissionPlan {
    std::vector<int> label;                        // Voronoi region per cell (-1 = outside search area)
    std::vector<Region> regions;
    std::vector<std::vector<Vec2>> sweepRoutes;    // One boustrophedon route per sweep drone
    std::vector<Vec2> acoTour;                     // Closed tour through the outer points of interest
};
