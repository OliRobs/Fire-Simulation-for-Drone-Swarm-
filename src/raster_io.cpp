#include "raster_io.h"
#include <gdal_priv.h>
#include <gdal_utils.h>
#include <ogr_spatialref.h>
#include <algorithm>
#include <cmath>
#include <initializer_list>
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
// North-up window in the local UTM zone
// ----------------------------------------------------------------------------
namespace {
    // Coordinate transformation between two WKT/SRS definitions, x/y in "GIS order" (lon, lat for degrees)
    struct Transform {
        OGRCoordinateTransformation* ct = nullptr;
        Transform(const OGRSpatialReference& from, const OGRSpatialReference& to)
            : ct(OGRCreateCoordinateTransformation(&from, &to)) {}
        ~Transform() { if (ct) OGRCoordinateTransformation::DestroyCT(ct); }
        Transform(const Transform&) = delete;
        Transform& operator=(const Transform&) = delete;
        bool apply(double& x, double& y) const { return ct && ct->Transform(1, &x, &y) != 0; }
    };
}

bool makeNorthUpWindow(const RasterGrid& full, double centreX, double centreY, double cellSize, RasterGrid& out,
                       int& utmZone, double& rotationDeg) {
    if (full.wkt.empty() || !full.projected) return false;

    OGRSpatialReference src, geo, utm;
    src.importFromWkt(full.wkt.c_str());
    geo.SetWellKnownGeogCS("NAD83");                       // LANDFIRE's datum
    for (OGRSpatialReference* r : { &src, &geo, &utm }) r->SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    // --- Pick the UTM zone from the centre's longitude (6-degree zones from 180 W) ---
    double lon = centreX, lat = centreY;
    if (!Transform(src, geo).apply(lon, lat)) return false;
    utmZone = std::clamp(static_cast<int>(std::floor((lon + 180.0) / 6.0)) + 1, 1, 60);
    utm.SetWellKnownGeogCS("NAD83");
    utm.SetUTM(utmZone, lat >= 0.0);
    utm.SetAxisMappingStrategy(OAMS_TRADITIONAL_GIS_ORDER);

    const Transform toUtm(src, utm), toSrc(utm, src);
    const double w = full.width * cellSize, h = full.height * cellSize;   // Window size on the ground

    // --- Place the window, moving it in from the file's edges if a corner falls outside ---
    // The north-up window is rotated relative to the file, so its corners are checked in the file's coordinates.
    double cx = 0, cy = 0;
    for (int attempt = 0; attempt < 10; ++attempt) {
        cx = centreX; cy = centreY;
        if (!toUtm.apply(cx, cy)) return false;
        cx = std::round(cx / cellSize) * cellSize;          // Keep cell edges on round coordinates
        cy = std::round(cy / cellSize) * cellSize;

        double lo[2] = { 1e300, 1e300 }, hi[2] = { -1e300, -1e300 };
        for (int corner = 0; corner < 4; ++corner) {
            double x = cx + ((corner & 1) ? w / 2 : -w / 2), y = cy + ((corner & 2) ? h / 2 : -h / 2);
            if (!toSrc.apply(x, y)) return false;
            lo[0] = std::min(lo[0], x); hi[0] = std::max(hi[0], x);
            lo[1] = std::min(lo[1], y); hi[1] = std::max(hi[1], y);
        }
        if (hi[0] - lo[0] > full.maxX - full.minX || hi[1] - lo[1] > full.maxY - full.minY) return false;   // Too small

        const double shiftX = std::max(full.minX - lo[0], 0.0) - std::max(hi[0] - full.maxX, 0.0);
        const double shiftY = std::max(full.minY - lo[1], 0.0) - std::max(hi[1] - full.maxY, 0.0);
        if (shiftX == 0.0 && shiftY == 0.0) break;           // Fits
        centreX += shiftX + (shiftX > 0 ? cellSize : shiftX < 0 ? -cellSize : 0.0);   // Small margin for rounding
        centreY += shiftY + (shiftY > 0 ? cellSize : shiftY < 0 ? -cellSize : 0.0);
    }

    // --- How far the file's grid is rotated from north here: map 1 km of UTM northing into the file ---
    double nx0 = cx, ny0 = cy, nx1 = cx, ny1 = cy + 1000.0;
    if (toSrc.apply(nx0, ny0) && toSrc.apply(nx1, ny1)) {
        rotationDeg = std::atan2(nx1 - nx0, ny1 - ny0) * 57.29577951308232;
    }

    char* wkt = nullptr;
    utm.exportToWkt(&wkt);
    out = full;
    out.wkt = wkt ? wkt : "";
    CPLFree(wkt);
    out.projected = true;
    out.minX = cx - w / 2; out.maxX = cx + w / 2;
    out.minY = cy - h / 2; out.maxY = cy + h / 2;
    out.originX = out.minX; out.originY = out.maxY;      // The new grid is its own pixel grid
    out.pixelW = out.pixelH = cellSize;
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
    args.AddString("-r");
    args.AddString(method == Resample::Mode ? "mode" : method == Resample::Max ? "max" : "average");

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
