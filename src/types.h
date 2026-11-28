#pragma once
#include <cmath>
#include <cstdint>
#include "config.h"

// Plain 2D vector in grid units (cell centres are at +0.5). Deliberately not sf::Vector2f, so the
// simulation and planning code builds and runs without SFML.
struct Vec2 {
    float x = 0.0f, y = 0.0f;
    constexpr Vec2() = default;
    constexpr Vec2(float x_, float y_) : x(x_), y(y_) {}
    constexpr Vec2 operator+(Vec2 o) const { return { x + o.x, y + o.y }; }
    constexpr Vec2 operator-(Vec2 o) const { return { x - o.x, y - o.y }; }
    constexpr Vec2 operator*(float s) const { return { x * s, y * s }; }
    constexpr Vec2 operator/(float s) const { return { x / s, y / s }; }
    Vec2& operator+=(Vec2 o) { x += o.x; y += o.y; return *this; }
};

inline float sqDist(Vec2 a, Vec2 b) {
    const float dx = a.x - b.x, dy = a.y - b.y;
    return dx * dx + dy * dy;
}

// Rectangle of grid cells [x0, x1) x [y0, y1) (half-open, so x1 - x0 is the width). Used to keep per-frame
// work to the part of the grid round the fire instead of all of it (GridSpec, grid.h, gives the grid's size).
struct GridRect {
    int x0 = 0, y0 = 0, x1 = 0, y1 = 0;   // Default: empty
    bool empty() const { return x1 <= x0 || y1 <= y0; }

    // Grows the rectangle to take in cell (x, y). An empty rectangle becomes that single cell.
    void include(int x, int y) {
        if (empty()) { x0 = x; y0 = y; x1 = x + 1; y1 = y + 1; return; }
        if (x < x0) x0 = x;
        if (y < y0) y0 = y;
        if (x >= x1) x1 = x + 1;
        if (y >= y1) y1 = y + 1;
    }

    // Smallest rectangle containing both this one and o (an empty one adds nothing)
    GridRect unite(const GridRect& o) const {
        if (empty()) return o;
        if (o.empty()) return *this;
        return { x0 < o.x0 ? x0 : o.x0, y0 < o.y0 ? y0 : o.y0, x1 > o.x1 ? x1 : o.x1, y1 > o.y1 ? y1 : o.y1 };
    }
};

enum class CellState : std::uint8_t { Unburned, Burning, Burned };

// One 10 m patch of landscape (SIM_CELL_SIZE_M): its terrain and LANDFIRE layers plus its fire state.
struct Cell {
    float elevation = 0.0f;  // Metres (LANDFIRE Elev)
    float slopeDeg = 0.0f;   // Slope steepness in degrees (LANDFIRE SlpD, or derived from elevation)
    float upslopeX = 0.0f;   // Unit vector pointing uphill, grid frame (x = east, y = south);
    float upslopeY = 0.0f;   //   (0, 0) on flat ground. Derived from elevation (the aspect).
    float burnTime = 0.0f;   // Minutes spent burning
    // Forest canopy (LANDFIRE CC, CH, CBH, CBD, converted to real units by the loader). All 0 = no trees.
    float canopyCoverPct = 0.0f;        // CC: percent of the ground covered by tree crowns
    float canopyHeightM = 0.0f;         // CH: height of the canopy top, m
    float canopyBaseHeightM = 0.0f;     // CBH: height of the canopy bottom (where fire can climb in), m
    float canopyBulkDensity = 0.0f;     // CBD: canopy fuel per unit volume, kg/m^3
    std::int16_t fuelModel = DEFAULT_FUEL_MODEL; // FBFM40 fuel model number (Scott & Burgan); 91-99 = nonburnable
    std::uint8_t roadClass = 0;                  // LANDFIRE road class 20-23, 0 = no road
    std::uint8_t fireType = 0;                   // How it burns if reached (FireType: 0 surface, 1 passive, 2 active crown); set by FireSimulation
    CellState state = CellState::Unburned;
};
