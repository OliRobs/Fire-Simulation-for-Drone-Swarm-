#include "raster_io.h"
#include <gdal_priv.h>
#include <gdal_utils.h>
#include <ogr_spatialref.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>

// ============================================================================
// Geospatial Raster I/O (GDAL)
//
// Every input layer (elevation, fuel, ...) is warped onto ONE simulation grid, normally defined by the DEM:
// same coordinate system, same extent, same cell size. That way a fuel map in a different projection or
// with a different footprint still lines up with the terrain, and a map that doesn't cover the terrain at
// all is rejected instead of being silently stretched over it.
// ============================================================================

// Registers all available GDAL drivers globally (must be called before loading files)
void registerRasterDrivers() { GDALAllRegister(); }

namespace {
    // RAII wrapper so datasets are always closed, even on early returns
    struct Dataset {
        GDALDatasetH h = nullptr;
        explicit Dataset(GDALDatasetH d) : h(d) {}
        ~Dataset() { if (h) GDALClose(h); }
        Dataset(const Dataset&) = delete;
        Dataset& operator=(const Dataset&) = delete;
    };

    // The dataset's coordinate system as WKT, or "" if it has none
    std::string wktOf(GDALDatasetH ds) {
        const OGRSpatialReference* srs = GDALDataset::FromHandle(ds)->GetSpatialRef();
        if (!srs) return "";
        char* wkt = nullptr;
        srs->exportToWkt(&wkt);
        std::string s = wkt ? wkt : "";
        CPLFree(wkt);
        return s;
    }

    // Full-precision text for a coordinate (std::to_string only keeps 6 decimals, too few for degrees)
    std::string fmt(double v) { return CPLSPrintf("%.17g", v); }
}

// ----------------------------------------------------------------------------
// Read a raster's footprint to use as the simulation grid
// ----------------------------------------------------------------------------
bool readRasterGrid(const char* filePath, int width, int height, RasterGrid& out) {
    Dataset ds(GDALOpen(filePath, GA_ReadOnly));
    if (!ds.h) {
        std::cerr << "[GDAL] Could not open " << filePath << "\n";
        return false;
    }

    double gt[6] = { 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 }; // Identity (pixel units) if the file has no geotransform
    if (GDALGetGeoTransform(ds.h, gt) != CE_None) {
        std::cerr << "[GDAL] " << filePath << " has no georeferencing; using its pixel extent\n";
    }
    if (gt[2] != 0.0 || gt[4] != 0.0) {
        std::cerr << "[GDAL] " << filePath << " is rotated; rotated rasters are not supported\n";
        return false;
    }

    // Geotransform: x = gt[0] + col * gt[1], y = gt[3] + row * gt[5] (gt[5] is negative for north-up)
    const int w = GDALGetRasterXSize(ds.h), h = GDALGetRasterYSize(ds.h);
    const double x0 = gt[0], x1 = gt[0] + w * gt[1];
    const double y0 = gt[3], y1 = gt[3] + h * gt[5];
    out.minX = std::min(x0, x1); out.maxX = std::max(x0, x1);
    out.minY = std::min(y0, y1); out.maxY = std::max(y0, y1);
    out.wkt = wktOf(ds.h);
    const OGRSpatialReference* srs = GDALDataset::FromHandle(ds.h)->GetSpatialRef();
    out.projected = srs && srs->IsProjected();
    out.width = width;
    out.height = height;
    out.originX = gt[0];
    out.originY = gt[3];
    out.pixelW = std::fabs(gt[1]);
    out.pixelH = std::fabs(gt[5]);

    return true;
}

// ----------------------------------------------------------------------------
// Cut a simulation window out of a large raster
// ----------------------------------------------------------------------------
bool makeWindow(const RasterGrid& full, double centreX, double centreY, double cellSize, RasterGrid& out) {
    const double w = full.width * cellSize, h = full.height * cellSize;   // Window size on the ground
    if (w > full.maxX - full.minX || h > full.maxY - full.minY) return false;

    // Keep the window inside the file: move the centre in from the edges if needed
    centreX = std::clamp(centreX, full.minX + w / 2, full.maxX - w / 2);
    centreY = std::clamp(centreY, full.minY + h / 2, full.maxY - h / 2);

    // Snap the left and top edges to the file's pixel boundaries, then clamp back inside after rounding
    double left = full.originX + std::round((centreX - w / 2 - full.originX) / full.pixelW) * full.pixelW;
    double top = full.originY - std::round((full.originY - (centreY + h / 2)) / full.pixelH) * full.pixelH;
    left = std::clamp(left, full.minX, full.maxX - w);
    top = std::clamp(top, full.minY + h, full.maxY);

    out = full;
    out.minX = left;
    out.maxX = left + w;
    out.maxY = top;
    out.minY = top - h;
    return true;
}

// ----------------------------------------------------------------------------
// Longitude/latitude -> grid coordinates
// ----------------------------------------------------------------------------
bool lonLatToGrid(const RasterGrid& grid, double lon, double lat, double& x, double& y) {
    if (grid.wkt.empty()) return false;
    OGRSpatialReference wgs84, target;
    wgs84.SetWellKnownGeogCS("WGS84");
    wgs84.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);   // x = longitude, y = latitude
    target.importFromWkt(grid.wkt.c_str());
    target.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    OGRCoordinateTransformation* ct = OGRCreateCoordinateTransformation(&wgs84, &target);
    if (!ct) return false;
    x = lon;
    y = lat;
    const bool ok = ct->Transform(1, &x, &y) != 0;
    OGRCoordinateTransformation::DestroyCT(ct);
    return ok;
}

// ----------------------------------------------------------------------------
// Warp one raster onto the simulation grid
// ----------------------------------------------------------------------------
bool loadRaster(const char* filePath, const RasterGrid& grid, Resample method, std::vector<float>& outData) {
    Dataset src(GDALOpen(filePath, GA_ReadOnly));
    if (!src.h) {
        std::cerr << "[GDAL] Could not open " << filePath << "\n";
        return false;
    }

    // Build gdalwarp arguments: in-memory output, Float32, NaN for "no data", target extent and size
    CPLStringList args;
    args.AddString("-of");  args.AddString("MEM");
    args.AddString("-ot");  args.AddString("Float32");
    args.AddString("-dstnodata"); args.AddString("nan");
    args.AddString("-te");
    args.AddString(fmt(grid.minX).c_str()); args.AddString(fmt(grid.minY).c_str());
    args.AddString(fmt(grid.maxX).c_str()); args.AddString(fmt(grid.maxY).c_str());
    args.AddString("-ts");
    args.AddString(std::to_string(grid.width).c_str()); args.AddString(std::to_string(grid.height).c_str());
    args.AddString("-r");   args.AddString(method == Resample::Mode ? "mode" : "average");

    // Reproject only when both sides have a coordinate system. If one is missing we can't check the
    // alignment, so say so and assume they already share one.
    const std::string srcWkt = wktOf(src.h);
    if (!grid.wkt.empty() && !srcWkt.empty()) {
        args.AddString("-t_srs"); args.AddString(grid.wkt.c_str());
    }
    else if (grid.wkt.empty() != srcWkt.empty()) {
        std::cerr << "[GDAL] Warning: " << filePath
                  << " and the terrain don't both have a coordinate system; assuming they match\n";
    }

    GDALWarpAppOptions* opts = GDALWarpAppOptionsNew(args.List(), nullptr);
    int usageError = 0;
    GDALDatasetH srcHandle = src.h;
    Dataset warped(GDALWarp("", nullptr, 1, &srcHandle, opts, &usageError));
    GDALWarpAppOptionsFree(opts);
    if (!warped.h) {
        std::cerr << "[GDAL] Could not warp " << filePath << " onto the terrain grid\n";
        return false;
    }

    // Read the warped band into the output buffer
    outData.assign(static_cast<size_t>(grid.width) * grid.height, std::numeric_limits<float>::quiet_NaN());
    GDALRasterBandH band = GDALGetRasterBand(warped.h, 1);
    if (GDALRasterIO(band, GF_Read, 0, 0, grid.width, grid.height, outData.data(),
                     grid.width, grid.height, GDT_Float32, 0, 0) != CE_None) {
        std::cerr << "[GDAL] Could not read " << filePath << "\n";
        return false;
    }

    // Reject files that don't cover the terrain at all (e.g. a US fuel map over a UK DEM)
    size_t valid = 0;
    for (float v : outData) valid += std::isnan(v) ? 0 : 1;
    if (valid == 0) {
        std::cerr << "[GDAL] " << filePath << " does not overlap the terrain extent; ignoring it\n";
        return false;
    }
    if (valid < outData.size()) {
        std::cerr << "[GDAL] Warning: " << filePath << " covers only "
                  << (100.0 * valid / outData.size()) << "% of the terrain\n";
    }

    std::cout << "[GDAL] Loaded " << filePath << "\n";
    return true;
}
