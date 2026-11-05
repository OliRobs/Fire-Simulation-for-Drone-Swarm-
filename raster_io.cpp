#include "raster_io.h"
#include <gdal_priv.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

// ============================================================================
// Geospatial Raster I/O Engine (GDAL Integration)
//
// Handles reading and resampling of GeoTIFF raster files (e.g., elevation, 
// fuel maps) into 1D simulation grid arrays using GDAL[cite: 1, 10, 12, 14].
// Includes automatic resampling to target grid dimensions and NoData/NaN 
// cleaning to ensure simulation stability[cite: 10].
// ============================================================================

// Registers all available GDAL drivers globally (must be called before loading files)[cite: 10]
void registerRasterDrivers() { GDALAllRegister(); }

// ----------------------------------------------------------------------------
// GeoTIFF Loader & Resampler
// ----------------------------------------------------------------------------
// Loads band 1 from a specified GeoTIFF raster file, resamples the spatial data 
// to match the simulation's (targetWidth x targetHeight) grid, and cleans 
// missing or invalid NoData sentinel values[cite: 10].
bool loadGeoTIFF(const char* filePath, std::vector<float>& outData,
    int targetWidth, int targetHeight, bool noDataToMin) {
    // Open the dataset in read-only mode[cite: 10]
    GDALDataset* ds = static_cast<GDALDataset*>(GDALOpen(filePath, GA_ReadOnly));
    if (!ds) {
        std::cerr << "[GDAL Error] Failed to open GeoTIFF file: " << filePath << "\n";
        return false;
    }

    // Access band 1 (primary data layer: elevation, fuel density, etc.)[cite: 10]
    GDALRasterBand* band = ds->GetRasterBand(1);
    int hasNoData = 0;
    // Extract NoData sentinel value (e.g., -9999) if defined in raster header[cite: 10]
    const float noData = static_cast<float>(band->GetNoDataValue(&hasNoData));

    // Resize output vector to hold resampled target grid cells[cite: 10]
    outData.resize(static_cast<size_t>(targetWidth) * targetHeight);

    // Configure GDAL RasterIO resampling arguments[cite: 10]
    GDALRasterIOExtraArg extra;
    INIT_RASTERIO_EXTRA_ARG(extra);
    extra.eResampleAlg = GRIORA_Average; // Bilinear/average resampling for continuous spatial rasters[cite: 10]

    // Read and resample raw GeoTIFF dimensions into output buffer array[cite: 10]
    CPLErr err = band->RasterIO(GF_Read, 0, 0, band->GetXSize(), band->GetYSize(),
        outData.data(), targetWidth, targetHeight,
        GDT_Float32, 0, 0, &extra);
    GDALClose(ds); // Close dataset handle[cite: 10]

    if (err != CE_None) {
        std::cerr << "[GDAL Error] Failed to read raster band data from: " << filePath << "\n";
        return false;
    }

    // Helper lambda to identify invalid/missing pixels (NaNs or header NoData sentinels)[cite: 10]
    auto invalid = [&](float v) { return std::isnan(v) || (hasNoData && v == noData); };

    // --- Clean and fill missing NoData values ---
    float minValid = std::numeric_limits<float>::max();
    for (float v : outData) {
        if (!invalid(v)) minValid = std::min(minValid, v);
    }

    // Choose replacement value: minimum valid value (useful for elevation rasters) or 0.0f[cite: 10]
    const float fill = (noDataToMin && minValid != std::numeric_limits<float>::max()) ? minValid : 0.0f;
    for (float& v : outData) {
        if (invalid(v)) v = fill;
    }

    std::cout << "[GDAL Success] Loaded raster map: " << filePath << "\n";
    return true;
}