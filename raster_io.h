#pragma once
#include "types.h"
#include <string>
#include <vector>

bool LoadLandfireDataset(const std::string& elev_path,
                         const std::string& slope_path,
                         const std::string& fuel_path,
                         const std::string& road_path,
                         std::vector<Cell>& grid,
                         int& width,
                         int& height);
