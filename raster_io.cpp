#include "raster_io.h"
#include <gdal_priv.h>
#include <iostream>

bool LoadLandfireDataset(const std::string& elev_path,
                         const std::string& slope_path,
                         const std::string& fuel_path,
                         const std::string& road_path,
                         std::vector<Cell>& grid,
                         int& width,
                         int& height) {
    GDALAllRegister();

    // 1. Open Elevation GeoTIFF to define primary dimensions
    GDALDataset* ds_elev = (GDALDataset*) GDALOpen(elev_path.c_str(), GA_ReadOnly);
    if (!ds_elev) {
        std::cerr << "Failed to open elevation file: " << elev_path << std::endl;
        return false;
    }

    width = ds_elev->GetRasterXSize();   // 8808
    height = ds_elev->GetRasterYSize();  // 6889
    const size_t total_pixels = static_cast<size_t>(width) * height;

    grid.clear();
    grid.resize(total_pixels);

    std::vector<float> buffer(total_pixels);

    // Read Elevation (meters)
    ds_elev->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, width, height, buffer.data(), width, height, GDT_Float32, 0, 0);
    for (size_t i = 0; i < total_pixels; ++i) {
        grid[i].elevation = buffer[i];
    }
    GDALClose(ds_elev);

    // 2. Read Slope (degrees)
    GDALDataset* ds_slope = (GDALDataset*) GDALOpen(slope_path.c_str(), GA_ReadOnly);
    if (ds_slope) {
        ds_slope->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, width, height, buffer.data(), width, height, GDT_Float32, 0, 0);
        for (size_t i = 0; i < total_pixels; ++i) {
            grid[i].slope_deg = buffer[i];
        }
        GDALClose(ds_slope);
    } else {
        std::cerr << "Warning: Could not open slope raster: " << slope_path << std::endl;
    }

    // 3. Read FBFM40 Fuel Models
    GDALDataset* ds_fuel = (GDALDataset*) GDALOpen(fuel_path.c_str(), GA_ReadOnly);
    if (ds_fuel) {
        ds_fuel->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, width, height, buffer.data(), width, height, GDT_Float32, 0, 0);
        for (size_t i = 0; i < total_pixels; ++i) {
            grid[i].fuel_model = static_cast<int>(buffer[i]);
        }
        GDALClose(ds_fuel);
    } else {
        std::cerr << "Warning: Could not open fuel model raster: " << fuel_path << std::endl;
    }

    // 4. Read Operational Roads
    GDALDataset* ds_road = (GDALDataset*) GDALOpen(road_path.c_str(), GA_ReadOnly);
    if (ds_road) {
        ds_road->GetRasterBand(1)->RasterIO(GF_Read, 0, 0, width, height, buffer.data(), width, height, GDT_Float32, 0, 0);
        for (size_t i = 0; i < total_pixels; ++i) {
            grid[i].is_road = (buffer[i] > 0.0f);
        }
        GDALClose(ds_road);
    } else {
        std::cerr << "Warning: Could not open operational roads raster: " << road_path << std::endl;
    }

    return true;
}
