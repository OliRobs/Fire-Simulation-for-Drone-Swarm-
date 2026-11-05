#pragma once
#include <vector>

// The only module that touches GDAL.
void registerRasterDrivers();

// Reads band 1, resamples (area-averaged) to targetWidth x targetHeight, and replaces NoData/NaN cells.
// noDataToMin: fill with the lowest valid value (for DEMs) instead of 0 (for fuel).
bool loadGeoTIFF(const char* filePath, std::vector<float>& outData,
                 int targetWidth, int targetHeight, bool noDataToMin);
