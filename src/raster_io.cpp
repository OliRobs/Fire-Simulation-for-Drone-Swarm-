#include "raster_io.h"
#include <gdal_priv.h>
#include <gdal_utils.h>
#include <ogr_spatialref.h>
#include <ogr_srs_api.h>
#include <cpl_conv.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
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

namespace {
    // ------------------------------------------------------------------------
    // Find PROJ's database (proj.db) and GDAL's data folder
    // ------------------------------------------------------------------------
    // Reprojecting the LANDFIRE layers (EPSG:5070 Albers) onto the DEM's grid (EPSG:26911 UTM) needs PROJ's
    // coordinate system database, proj.db. PROJ looks for it in the folder named by the PROJ_DATA (PROJ 9) or
    // PROJ_LIB (older) environment variable, and a Windows GDAL unpacked into a plain folder sets neither. Without
    // it every reprojection fails ("Cannot find proj.db", "Cannot find coordinate operations"), so the fuel,
    // roads and canopy layers are all dropped. If the variables aren't set, look in the usual places instead:
    // the GISInternals layout used for this project (C:\gdal\bin\proj9\share, as in UAVFireMap.vcxproj), a
    // conda-style C:\gdal\share\proj, next to the program, and the Linux system folders.
    // ------------------------------------------------------------------------
    // Value of an environment variable, or "" if it isn't set.
    // MSVC flags std::getenv as unsafe (warning C4996, an error with the project's SDL checks on), so on MSVC
    // this uses _dupenv_s, which returns a copy the caller must free; elsewhere plain std::getenv.
    // ------------------------------------------------------------------------
    std::string envVar(const char* name) {
#ifdef _MSC_VER
        char* value = nullptr;
        size_t len = 0;
        std::string out;
        if (_dupenv_s(&value, &len, name) == 0 && value) out = value;
        std::free(value);
        return out;
#else
        const char* value = std::getenv(name);
        return value ? value : "";
#endif
    }

    std::string firstFolderWith(std::initializer_list<const char*> folders, const char* file) {
        for (const char* f : folders) {
            std::error_code ec;
            if (f && *f && std::filesystem::exists(std::filesystem::path(f) / file, ec)) return f;
        }
        return {};
    }

    void configureProjAndGdalData() {
        // --- PROJ database ---
        const std::string projData = envVar("PROJ_DATA"), projLib = envVar("PROJ_LIB");
        std::string projDir = firstFolderWith({ projData.c_str(), projLib.c_str() }, "proj.db");
        if (projDir.empty()) {
            projDir = firstFolderWith({ "C:/gdal/bin/proj9/share", "C:/gdal/bin/proj8/share", "C:/gdal/bin/proj7/share",
                                        "C:/gdal/share/proj", "proj9/share", "proj/share", "share/proj",
                                        "/usr/share/proj", "/usr/local/share/proj" }, "proj.db");
            if (!projDir.empty()) {
                const char* const paths[] = { projDir.c_str(), nullptr };
                OSRSetPROJSearchPaths(paths);   // Same effect as setting PROJ_DATA, for this program only
                std::cout << "[GDAL] PROJ database (proj.db) found in " << projDir << "\n";
            }
            else {
                std::cerr << "[GDAL] WARNING: PROJ database (proj.db) not found, so the LANDFIRE layers can't be reprojected\n"
                             "       onto the terrain grid and will be skipped. Set the PROJ_DATA environment variable to the\n"
                             "       folder containing proj.db (e.g. C:\\gdal\\bin\\proj9\\share).\n";
            }
        }

        // --- GDAL support files (CSV tables, etc.): GDAL_DATA ---
        if (!CPLGetConfigOption("GDAL_DATA", nullptr)) {
            const std::string gdalData = firstFolderWith({ "C:/gdal/bin/gdal-data", "C:/gdal/share/gdal", "gdal-data",
                                                           "/usr/share/gdal", "/usr/local/share/gdal" }, "gdalvrt.xsd");
            if (!gdalData.empty()) CPLSetConfigOption("GDAL_DATA", gdalData.c_str());
        }
    }
}

// Registers all available GDAL drivers globally (must be called before loading files), after pointing PROJ and
// GDAL at their data folders (see configureProjAndGdalData)
void registerRasterDrivers() {
    configureProjAndGdalData();
    GDALAllRegister();
}

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
    // USGS DEMs carry a compound coordinate system (UTM + NAVD88 heights). Keep only the horizontal part:
    // the grid is a 2D layout, and a vertical CRS in the warp target can trigger height conversions.
    if (srs && srs->IsCompound()) {
        OGRSpatialReference horizontal(*srs);
        horizontal.StripVertical();
        char* wkt = nullptr;
        horizontal.exportToWkt(&wkt);
        out.wkt = wkt ? wkt : "";
        CPLFree(wkt);
    }
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

// ----------------------------------------------------------------------------
// Mosaic of DEM tiles
// ----------------------------------------------------------------------------
bool buildMosaic(const std::vector<std::string>& tiles, const std::string& vrtPath) {
    if (tiles.empty()) return false;
    std::vector<const char*> names;
    for (const std::string& t : tiles) names.push_back(t.c_str());

    // gdalbuildvrt equivalent. Where tiles overlap (USGS tiles have a few metres of buffer) the values are
    // the same, so it doesn't matter which tile is read.
    int usageError = 0;
    Dataset vrt(GDALBuildVRT(vrtPath.c_str(), static_cast<int>(names.size()), nullptr, names.data(), nullptr,
                             &usageError));
    if (!vrt.h) {
        std::cerr << "[GDAL] Could not join the DEM tiles into one mosaic (different coordinate systems?)\n";
        return false;
    }
    std::cout << "[GDAL] DEM mosaic of " << tiles.size() << " tile(s): " << GDALGetRasterXSize(vrt.h) << " x "
              << GDALGetRasterYSize(vrt.h) << " pixels\n";
    return true;   // Closing the handle writes the VRT to vrtPath, ready to be opened again
}

// ----------------------------------------------------------------------------
// Grid-cell terrain from the full-resolution DEM
// ----------------------------------------------------------------------------
bool loadFineTopography(const char* demPath, const RasterGrid& grid, std::vector<CellTopography>& out,
                        double& finePixelM) {
    Dataset ds(GDALOpen(demPath, GA_ReadOnly));
    if (!ds.h) {
        std::cerr << "[GDAL] Could not open " << demPath << "\n";
        return false;
    }
    double gt[6];
    if (GDALGetGeoTransform(ds.h, gt) != CE_None || gt[2] != 0.0 || gt[4] != 0.0) {
        std::cerr << "[GDAL] " << demPath << " has no north-up georeferencing\n";
        return false;
    }
    GDALRasterBandH band = GDALGetRasterBand(ds.h, 1);
    const int rasterW = GDALGetRasterXSize(ds.h), rasterH = GDALGetRasterYSize(ds.h);
    int hasNoData = 0;
    const double noData = GDALGetRasterNoDataValue(band, &hasNoData);

    // --- Fine pixel size: the DEM's, nudged so a whole number of pixels fits in one cell ---
    // (1 m DEM, 10 m cells: exactly 10 x 10 pixels per cell, read 1:1 with no resampling)
    const double cell = grid.cellSizeX();
    const double demPixel = std::fabs(gt[1]);
    const int factor = std::max(1, static_cast<int>(std::lround(cell / demPixel)));
    finePixelM = cell / factor;

    // Source pixel coordinates of the buffer's top-left corner (the buffer has a one-pixel halo)
    const double colStart = (grid.minX - finePixelM - gt[0]) / gt[1];
    const double scale = finePixelM / demPixel;          // Source pixels per buffer pixel (1.0 for 1 m)
    const bool oneToOne = std::fabs(scale - 1.0) < 1e-6
        && std::fabs(colStart - std::round(colStart)) < 1e-3
        && std::fabs((gt[3] - grid.maxY - finePixelM) / std::fabs(gt[5])
                     - std::round((gt[3] - grid.maxY - finePixelM) / std::fabs(gt[5]))) < 1e-3;

    // Reads buffer pixels [i0, i1) x [j0, j1) of the strip whose first buffer row starts at source row
    // rowStart, into the right place in `buf` (width bufW). Parts outside the raster stay NaN.
    auto readPart = [&](double rowStart, int bufW, int bufH, std::vector<float>& buf) -> bool {
        // Which buffer pixels fall inside the raster (each buffer pixel spans `scale` source pixels)
        const int i0 = std::max(0, static_cast<int>(std::ceil(-colStart / scale - 1e-9)));
        const int i1 = std::min(bufW, static_cast<int>(std::floor((rasterW - colStart) / scale + 1e-9)));
        const int j0 = std::max(0, static_cast<int>(std::ceil(-rowStart / scale - 1e-9)));
        const int j1 = std::min(bufH, static_cast<int>(std::floor((rasterH - rowStart) / scale + 1e-9)));
        if (i1 <= i0 || j1 <= j0) return true;           // Strip entirely outside the DEM
        const int partW = i1 - i0, partH = j1 - j0;
        std::vector<float> part(static_cast<size_t>(partW) * partH);

        // Source window for that part. 1:1 reads are plain copies; otherwise GDAL resamples bilinearly from a
        // floating-point window, so cells still line up with the ground exactly.
        const double dfX = colStart + i0 * scale, dfY = rowStart + j0 * scale;
        const double dfW = partW * scale, dfH = partH * scale;
        GDALRasterIOExtraArg arg;
        INIT_RASTERIO_EXTRA_ARG(arg);
        int nX = static_cast<int>(std::lround(dfX)), nY = static_cast<int>(std::lround(dfY));
        int nW = partW, nH = partH;
        if (!oneToOne) {
            arg.eResampleAlg = GRIORA_Bilinear;
            arg.bFloatingPointWindowValidity = TRUE;
            arg.dfXOff = dfX; arg.dfYOff = dfY; arg.dfXSize = dfW; arg.dfYSize = dfH;
            nX = static_cast<int>(std::floor(dfX));
            nY = static_cast<int>(std::floor(dfY));
            nW = std::min(rasterW, static_cast<int>(std::ceil(dfX + dfW - 1e-9))) - nX;
            nH = std::min(rasterH, static_cast<int>(std::ceil(dfY + dfH - 1e-9))) - nY;
        }
        nX = std::clamp(nX, 0, rasterW - 1);
        nY = std::clamp(nY, 0, rasterH - 1);
        nW = std::clamp(nW, 1, rasterW - nX);
        nH = std::clamp(nH, 1, rasterH - nY);
        if (GDALRasterIOEx(band, GF_Read, nX, nY, nW, nH, part.data(), partW, partH, GDT_Float32, 0, 0, &arg)
            != CE_None) {
            return false;
        }
        // Copy into the strip buffer, turning NoData into NaN
        for (int j = 0; j < partH; ++j) {
            for (int i = 0; i < partW; ++i) {
                float v = part[static_cast<size_t>(j) * partW + i];
                if ((hasNoData && v == static_cast<float>(noData)) || v <= -9000.0f) v = NAN;
                buf[static_cast<size_t>(j0 + j) * bufW + (i0 + i)] = v;
            }
        }
        return true;
    };

    // --- Read strip by strip and summarise each strip's cells ---
    constexpr int STRIP_CELL_ROWS = 50;                  // 50 cell rows: (1200 x 10 + 2) x 502 floats = 24 MB
    const int bufW = grid.width * factor + 2;
    out.assign(static_cast<size_t>(grid.width) * grid.height, CellTopography{});
    std::vector<float> buf;
    std::vector<CellTopography> stripOut;
    int lastReported = -1;
    for (int r0 = 0; r0 < grid.height; r0 += STRIP_CELL_ROWS) {
        const int rows = std::min(STRIP_CELL_ROWS, grid.height - r0);
        const int bufH = rows * factor + 2;
        buf.assign(static_cast<size_t>(bufW) * bufH, NAN);
        // Source row of the strip's first buffer row (the halo row just above cell row r0)
        const double rowStart = (gt[3] - (grid.maxY - r0 * cell + finePixelM)) / std::fabs(gt[5]);
        if (!readPart(rowStart, bufW, bufH, buf)) {
            std::cerr << "[GDAL] Could not read " << demPath << "\n";
            return false;
        }
        aggregateFineTopography(buf, grid.width, rows, factor, static_cast<float>(finePixelM),
                                static_cast<float>(finePixelM), stripOut);
        std::copy(stripOut.begin(), stripOut.end(), out.begin() + static_cast<size_t>(r0) * grid.width);

        // Progress: the window is 120 million 1 m pixels, so this takes a little while
        const int pct = 100 * (r0 + rows) / grid.height;
        if (pct / 10 != lastReported) {
            lastReported = pct / 10;
            std::cout << "\r[GDAL] Reading the DEM at " << finePixelM << " m: " << pct << "%" << std::flush;
        }
    }
    std::cout << "\n";

    size_t covered = 0;
    for (const CellTopography& t : out) covered += t.validPixels > 0 ? 1 : 0;
    if (covered == 0) {
        std::cerr << "[GDAL] " << demPath << " does not overlap the simulation grid\n";
        return false;
    }
    std::cout << "[GDAL] DEM covers " << (100.0 * covered / out.size()) << "% of the grid, " << factor << " x "
              << factor << " pixels per cell\n";
    return true;
}
