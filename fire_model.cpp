#include "fire_model.h"
#include <cmath>
#include <algorithm>

FireSimulation::FireSimulation(std::vector<Cell>& grid, int width, int height)
    : grid_(grid), width_(width), height_(height) {}

bool FireSimulation::IsNonBurnable(int fuel_code, bool is_road) const {
    if (is_road) return true; // Operational roads serve as non-burnable barriers
    // Non-burnable LANDFIRE FBFM40 codes: 91 (Urban), 92 (Snow/Ice), 93 (Ag), 98 (Water), 99 (Barren)
    return (fuel_code == 91 || fuel_code == 92 || fuel_code == 93 || fuel_code == 98 || fuel_code == 99 || fuel_code <= 0);
}

float FireSimulation::CalculateCellROS(const Cell& cell, float wind_speed, float wind_dir) const {
    if (IsNonBurnable(cell.fuel_model, cell.is_road)) {
        return 0.0f;
    }

    // Convert slope degrees to percent slope: percent_slope = tan(slope_deg) * 100
    float slope_rad = cell.slope_deg * (3.14159265f / 180.0f);
    float slope_percent = std::tan(slope_rad) * 100.0f;

    // Slope influence factor phi_s
    float phi_s = 5.275f * std::pow(slope_percent / 100.0f, 2.0f);

    // Base Rate of Spread (m/s) lookup by Scott & Burgan categories
    float base_ros = 0.2f;
    if (cell.fuel_model >= 101 && cell.fuel_model <= 109)      base_ros = 0.6f; // Low/High Grass (GR)
    else if (cell.fuel_model >= 121 && cell.fuel_model <= 124) base_ros = 0.4f; // Low/High Shrub (SH)
    else if (cell.fuel_model >= 141 && cell.fuel_model <= 149) base_ros = 0.5f; // Timber-Shrub (TS)
    else if (cell.fuel_model >= 161 && cell.fuel_model <= 165) base_ros = 0.3f; // Timber-Understory (TU)
    else if (cell.fuel_model >= 181 && cell.fuel_model <= 189) base_ros = 0.2f; // Timber Litter (TL)

    // Combined Rate of Spread
    return base_ros * (1.0f + phi_s + (0.1f * wind_speed));
}

void FireSimulation::Ignite(int x, int y) {
    int idx = y * width_ + x;
    if (idx >= 0 && idx < static_cast<int>(grid_.size())) {
        if (!IsNonBurnable(grid_[idx].fuel_model, grid_[idx].is_road)) {
            grid_[idx].state = 1; // Burning
        }
    }
}

void FireSimulation::Update(float delta_time, float wind_speed, float wind_dir) {
    std::vector<int> to_ignite;

    for (int y = 0; y < height_; ++y) {
        for (int x = 0; x < width_; ++x) {
            int idx = y * width_ + x;
            if (grid_[idx].state == 1) { // Burning cell
                grid_[idx].burn_time += delta_time;

                // Simple 4-neighbor spread propagation
                const int dx[] = {0, 0, -1, 1};
                const int dy[] = {-1, 1, 0, 0};

                for (int i = 0; i < 4; ++i) {
                    int nx = x + dx[i];
                    int ny = y + dy[i];

                    if (nx >= 0 && nx < width_ && ny >= 0 && ny < height_) {
                        int nidx = ny * width_ + nx;
                        if (grid_[nidx].state == 0) {
                            float n_ros = CalculateCellROS(grid_[nidx], wind_speed, wind_dir);
                            if (n_ros > 0.0f && grid_[idx].burn_time > (30.0f / n_ros)) {
                                to_ignite.push_back(nidx);
                            }
                        }
                    }
                }

                if (grid_[idx].burn_time > 300.0f) {
                    grid_[idx].state = 2; // Burned out
                }
            }
        }
    }

    for (int nidx : to_ignite) {
        grid_[nidx].state = 1;
    }
}
