#pragma once
#include "types.h"
#include <vector>

class FireSimulation {
public:
    FireSimulation(std::vector<Cell>& grid, int width, int height);

    void Update(float delta_time, float wind_speed, float wind_dir);
    void Ignite(int x, int y);

private:
    std::vector<Cell>& grid_;
    int width_;
    int height_;

    bool IsNonBurnable(int fuel_code, bool is_road) const;
    float CalculateCellROS(const Cell& cell, float wind_speed, float wind_dir) const;
};
