#pragma once
#include <vector>
#include "types.h"

// Fills each cell's upslope direction (the aspect, needed to add slope and wind as vectors in the Rothermel
// model) from the elevation grid. If haveSlopeLayer is false, slopeDeg is derived from elevation too;
// otherwise the measured LANDFIRE slope is kept and only the direction comes from the DEM.
// cellSizeX/Y: ground size of one cell in metres, the same units as elevation.
void computeTopography(std::vector<Cell>& grid, float cellSizeX, float cellSizeY, bool haveSlopeLayer);

// Terrain summary of one simulation cell, worked out from the finer DEM pixels inside it
struct CellTopography {
    float elevation = 0.0f;   // Mean elevation of the pixels, m
    float slopeDeg = 0.0f;    // Mean of the pixels' own slopes, degrees (keeps small steep features)
    float upslopeX = 0.0f;    // Unit vector uphill (x = east, y = south) from the mean gradient; (0, 0) if flat
    float upslopeY = 0.0f;
    int validPixels = 0;      // Pixels with data; 0 = the DEM doesn't cover this cell
};

// Summarises a block of fine DEM pixels into simulation cells, factor x factor pixels per cell.
// Used for the USGS 1 m DEM: slope is worked out on every 1 m pixel (Horn's method, as above) and then
// averaged into the 10 m cell, so a cell crossing a gully or a steep bank keeps that steepness. Slope taken
// from 10 m averaged elevations would smooth those features away.
//   fine:      (cellsX * factor + 2) x (cellsY * factor + 2) elevations, row-major, north at the top, with a
//              ONE-PIXEL HALO all round (the pixels just outside the block, so edge pixels get a full 3x3).
//              NaN = no data. A NaN neighbour is replaced by the centre pixel's own value.
//   finePixelX/Y: ground size of one fine pixel, m.
//   out:       cellsX x cellsY results, row-major.
void aggregateFineTopography(const std::vector<float>& fine, int cellsX, int cellsY, int factor,
                             float finePixelX, float finePixelY, std::vector<CellTopography>& out);
