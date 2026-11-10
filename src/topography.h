#pragma once
#include <vector>
#include "types.h"

// Fills each cell's upslope direction (the aspect, needed to add slope and wind as vectors in the Rothermel
// model) from the elevation grid. If haveSlopeLayer is false, slopeDeg is derived from elevation too;
// otherwise the measured LANDFIRE slope is kept and only the direction comes from the DEM.
// cellSizeX/Y: ground size of one cell in metres, the same units as elevation.
void computeTopography(std::vector<Cell>& grid, float cellSizeX, float cellSizeY, bool haveSlopeLayer);
